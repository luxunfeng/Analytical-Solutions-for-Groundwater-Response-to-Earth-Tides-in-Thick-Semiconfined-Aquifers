%% For tidal behavior of a well in relatively leaky reservoirs
% keep k_overburden/k/b is constant
% change the values of C_D and S, and change hD to see the 5% variation
% surface. it is noted that C is linear with h for vertical well.
% subroutine needed: solve_tan_eq_n.m
% vertical well
clc;close all;clear
h = [20];
fai=0.25;
ct=1.02e-9;
rw=0.1;
xw = 0; x = rw;
yw = 0; y = 0;

z = 0; % the location of measurement gauge

kr=4e-12;kx=kr;ky=kr;kz=kr;
mu=0.003;
% isotropic permeability
etax=kx/(fai*mu*ct);etaz=kz/(fai*mu*ct);etay=ky/(fai*mu*ct);
% typical wellbore storage for horizontal well
C = 0.0001*h; %bbl/psi, C = cwVw, depends on wellbore geometry and
% compressibility of water. should not be a variable.
C = C / 43366.7; % m3/Pa
Cdd = C./(2*pi*h*fai*ct*rw^2); %non-dimensional wellbore storage

% nondimensional form
rd = rw / rw;
hdd=h/(rw);
zd = z/(rw);

period_M2 = 12.421*3600;
period_K1 = 23.934*3600;
period_Mf = 13.661*24*3600;
period_Mm = 27.555*24*3600;
period_Ssa = 0.5*365*24*3600;

period_components = [period_M2]; %s
nondim_co = kr/(fai*mu*ct*(rw)^2);
period_components_nondim = nondim_co * period_components;
omega_components_non = 2*pi./period_components_nondim * sqrt(-1);


S = [-2:0.2:5];

%%%%%%%%%%%%%% semi confined reservoir
color = {'-r','-b','-g','-k'};
color2 = {'--r','--b','--g','--k'};

%% overburden layer
b = 100;
bd = b / rw;
k_overburden = 4e-13;
Hdd = k_overburden / kr / bd ./ hdd;
criterion = hdd.*sqrt(Hdd);

%% previous model
fig = figure(1);
for ii = 1 : length(S)
    Hd = Hdd;
    hd = hdd;
    C_D = Cdd;

    omega = omega_components_non /sqrt(-1) * nondim_co;

    Sd = 1/ (2 * C_D);
    Td = 2 * pi ^2 * kr * h ./(C * mu * omega);
    ad = sqrt(2 * pi * sqrt(-1) * Sd/Td);
    betad = sqrt(Hd + ad^2);

    H = (ad/betad)^2 / (1 + ad^2/ (2 * Sd * betad) * besselk(0, betad) / besselk(1, betad) ...
        + pi * sqrt(-1) * S(ii) /Td);

    amplitude1(ii) = double(abs(H));
    phase1(ii) =  double(angle(H))/pi * 180;
    value1(ii) = pi * sqrt(-1) * S(ii) /Td;
    value2(ii) = ad^2/ (2 * Sd * betad) * besselk(0, betad) / besselk(1, betad);


end

yyaxis right
plot(S, amplitude1)

yyaxis left
plot(S, phase1)


figure(2)
plot(S, abs(value1),'-r'); hold on
plot(S, abs(value2),'-b')


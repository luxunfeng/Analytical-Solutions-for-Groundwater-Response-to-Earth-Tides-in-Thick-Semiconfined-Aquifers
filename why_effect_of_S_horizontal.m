%% For tidal behavior of a well in relatively leaky reservoirs
% keep k_overburden/k/b is constant
% change the values of C_D and S, and change hD to see the 5% variation
% surface. it is noted that C is linear with L for horizontal well.
% if well length constant, then C_D constant.
% here, keep L constant and change h, so hD changes.
% maybe could keep h constant and change L, so hD and CD change
% subroutine needed: solve_tan_eq_n.m
% Horizontal well
clc;close all;clear
%% parameters
h = 30;
fai=0.25;
ct=1.02e-9;
l=300; %length of horizontal well
xf=l;
rw=0.1;
xw=0;x=0;
zw=0.5* h;z=zw+rw; %look at the property of centering location of well bore（0,0,zw+rw）
yw=0;y=0;

kr=4e-12;kx=kr;ky=kr;kz=kr;
mu=0.003;
% isotropic permeability
etax=kx/(fai*mu*ct);etaz=kz/(fai*mu*ct);etay=ky/(fai*mu*ct);
% typical wellbore storage for horizontal well
C = 3.3333e-04*l; %bbl/psi
C = C / 43366.7; % m3/Pa
C_D = C/(2*pi*l*fai*ct*(l/2)^2); %non-dimensional wellbore storage

%% nondimensional form
hdd=h/(l/2);
ld=l/(l/2);
zdd = z/(l/2);
zwdd = zw/(l/2);
rwd = rw/(l/2);

period_M2 = 12.421*3600;
period_K1 = 23.934*3600;
period_Mf = 13.661*24*3600;
period_Mm = 27.555*24*3600;
period_Ssa = 0.5*365*24*3600;

period_components = [period_M2]; %s
nondim_co = kr/(fai*mu*ct*(l/2)^2);
period_components_nondim = nondim_co * period_components;
omega_components_non = 2*pi./period_components_nondim*sqrt(-1);

S = [-2:0.2:5]; 

%%%%%%%%%%%%%% semi confined reservoir
color = {'-r','-b','-g','-k'};
color2 = {'--r','--b','--g','--k'};

%% overburden layer
b = 100;
bd = b/(l/2);
k_overburden = 4e-13;
Hdd = k_overburden / kr / bd ./ hdd;

criterion = hdd.*sqrt(Hdd);

fig = figure(1);
%% previous model
for ii = 1 : length(S)
    Hd = Hdd;
    hd = hdd;
    zd = zdd;
    zwd = zwdd;
    s = omega_components_non;
    sum = 0;
    for j = 1 : 120
        sum = sum + asinh(pi/2/sqrt(s + Hd + (j*pi/hd)^2))*cos(j*pi*zd/hd)*cos(j*pi*zwd/hd);
    end
    A_bar =4/s * 1/ hd/2 *asinh(pi/2/sqrt(s + Hd)) +  8/s * 1 / hd / 2 * sum;
    A = s /(s + Hd);
    H = A/(1 + C_D * S(ii) * s + C_D * s^2 * A_bar);

    imag1(ii) = imag(C_D * S(ii) * s);
    imag2(ii) = imag(C_D * s^2 * A_bar);
    amplitude(ii) = abs(H);
    phase(ii) = angle(H);

end
   
plot(S, imag1);
hold on
plot(S, imag2)

figure(2)
yyaxis left
plot(S, amplitude);
yyaxis right
plot(S, phase);
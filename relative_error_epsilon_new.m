%% For tidal behavior of a well in relatively leaky reservoirs
% keep k_overburden is constant
% change the values of h, zw, z
% subroutine needed: solve_tan_eq_n.m
% vertical well
clc;close all;clear
% bh = 6000;
h = 10:1:1000;
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
C = 0.0001*h; %bbl/psi

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

legend_list = {'M2','K1','Mf','Mm','Ssa'}

period_components = [period_M2, period_K1, period_Mf, period_Mm, period_Ssa]; %s
nondim_co = kr/(fai*mu*ct*(rw)^2);
period_components_nondim = nondim_co * period_components;
omega_components_non = 2*pi./period_components_nondim *sqrt(-1);


S = 1;

%%%%%%%%%%%%%% semi confined reservoir
color = {'-r','-b','-g','-k'};
color2 = {'--r','--b','--g','--k'};
linewidth = 2;

%% overburden layer
% b = bh./h;
b=400;
bd = b / rw;
k_overburden = 4e-13;
Hdd = k_overburden / kr ./ bd ./ hdd;

%% previous model
for jj = 1 : length(Hdd)
    Hd = Hdd(jj);
    hd = hdd(jj);
    C_D = Cdd(jj);

    for i = 1 : length(omega_components_non)
        omega = omega_components_non(i) /sqrt(-1) * nondim_co;

        Sd = 1/ (2 * C_D);
        Td = 2 * pi ^2 * kr * h(jj) ./(C(jj) * mu * omega);
        ad = sqrt(2 * pi * sqrt(-1) * Sd/Td);
        betad = sqrt(Hd + ad^2);

        H = (ad/betad)^2 / (1 + ad^2/ (2 * Sd * betad) * besselk(0, betad) / besselk(1, betad) ...
            + pi * sqrt(-1) * S /Td);

        amplitude1(jj,i) = double(abs(H));
        phase1(jj,i) =  double(angle(H))/pi * 180;
    end
    %% new model
    lambda_d = Hd * hd;

    x = solve_tan_eq_n(hd*lambda_d,200);
    xx = x/hd;

    for i = 1 : length(omega_components_non)

        s = omega_components_non(i);

        sum1 = 0;
        for j = 1 : 200
            sum1 = sum1 + 1/xx(j) * sin(xx(j) * hd) * cos(xx(j)*zd)...
                * besselk(0, sqrt(s + xx(j)^2))...
                * (xx(j)^2 + lambda_d^2)/(hd * (xx(j)^2 + lambda_d^2) + lambda_d);
        end
        sum1 = 2 * sum1;

        sum2 = 0;
        for j = 1 : 200
            sum2 = sum2 - lambda_d*2/ (s + xx(j)^2)*cos(xx(j)*zd)*cos(xx(j)*hd)...
                * (xx(j)^2 + lambda_d^2)/(hd * (xx(j)^2 + lambda_d^2) + lambda_d);
        end

        H = (sum2 + 1)/(1 + C_D * S * s + C_D * s * sum1);
        amplitude2(jj,i) = double(abs(H));
        phase2(jj,i) =  double(angle(H))/pi * 180;
    end
end

biot = hdd.*sqrt(Hdd);
color = {'-r','-b','-g','-k','-m'};
color2 = {'--r','--b','--g','--k','--m'};
linewidth = 2;
subplot(1,2,1)
for ii = 1 : length(omega_components_non)

    plot(biot,100*abs(amplitude1(:,ii)-amplitude2(:,ii))./amplitude2(:,ii),color{ii},'LineWidth',linewidth)
    hold on
    xlabel('$B_{LSR}$','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|\frac{|H_1|-|H_2|}{|H_2|}| \%$','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    xlim([0,0.3])
    grid on
    % title2 = ["$S = 1, C_D ="+ string(C_D)+"$"];
    % title(title2,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
    ax = gca;
    set(ax, 'FontSize', 20);
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
end


subplot(1,2,2)
for ii = 1 : length(omega_components_non)

    plot(biot,100*abs(phase1(:,ii) - phase2(:,ii))./abs(phase2(:,ii)),color{ii},'LineWidth',linewidth)
    hold on
    xlabel('$B_{LSR}$','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}| \% $','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
    xlim([0,0.3])
    ax = gca;
    set(ax, 'FontSize', 20);
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    legend(legend_list)
end





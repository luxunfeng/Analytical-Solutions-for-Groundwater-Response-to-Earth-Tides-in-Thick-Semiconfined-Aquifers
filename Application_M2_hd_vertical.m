%% For tidal behavior of a well in relatively leaky reservoirs
% k_overburden changes
% change the values of h, zw, z
% subroutine needed: solve_tan_eq_n.m
% vertical well
clc;close all;clear
h = 48;
rho = 1000;
g = 10;
SS = 2.7E-5;
T = 9.6e-6;
fai_ct = SS/(h*rho*g)


rw=0.11;
rc = 0.0365;
xw = 0; x = rw;
yw = 0; y = 0;

z = 0; % the location of measurement gauge
mu=0.001;
kr=T/h*mu/rho/g;
kx=kr;ky=kr;kz=kr;

% isotropic permeability
etax=kx/(fai_ct*mu);etaz=kz/(fai_ct*mu);etay=ky/(fai_ct*mu);
% typical wellbore storage for horizontal well
C = pi * rc^2/ rho / g; % m3/Pa
Cdd = C./(2*pi*h*fai_ct*rw^2); %non-dimensional wellbore storage

% nondimensional form
rd = rw / rw;
hdd=h/(rw);
zd = z/(rw);

period_M2 = 12.421*3600;


period_components = [period_M2]; %s
nondim_co = kr/(fai_ct*mu*(rw)^2);
period_components_nondim = nondim_co * period_components;
omega_components_non = 2*pi./period_components_nondim*sqrt(-1);


S = 0; 

%%%%%%%%%%%%%% semi confined reservoir
color = {'-r','-b','-g','-k'};
color2 = {'--r','--b','--g','--k'};

%% overburden layer
b = 277;
bd = b / rw;
k_overburden = [0.01:0.01:2]*kr;
Hdd = k_overburden / kr / bd ./ hdd;

%% previous model
for jj = 1 : length(Hdd)
    Hd = Hdd(jj);
    hd = hdd(1);
    C_D = Cdd(1);

    
        omega = omega_components_non(1) /sqrt(-1) * nondim_co;

        Sd = 1/ (2 * C_D);
        Td = 2 * pi ^2 * kr * h ./(C * mu * omega);
        ad = sqrt(2 * pi * sqrt(-1) * Sd/Td);
        betad = sqrt(Hd + ad^2);

        H = (ad/betad)^2 / (1 + ad^2/ (2 * Sd * betad) * besselk(0, betad) / besselk(1, betad) ...
            + pi * sqrt(-1) * S /Td);

        amplitude1(jj) = double(abs(H));
        phase1(jj) =  double(angle(H))/pi * 180;
    
%% new model
    lambda_d = Hd * hd;

    x = solve_tan_eq_n(hd*lambda_d,200);
    xx = x/hd;

  
        
        s = omega_components_non(1);

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
        amplitude2(jj) = double(abs(H));
        phase2(jj) =  double(angle(H))/pi * 180;
 

    
  
end

  fig = figure(1);
    fig.Position = [100 100 800 600];
    subplot(2,2,1)
    semilogx(Hdd,amplitude1,color{1},'LineWidth',1)
    hold on
    semilogx(Hdd,amplitude2,color2{1},'LineWidth',1)
    xlabel('$H_{D}$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|H(i\omega_{D:M2})|$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
  
    grid on

    title1 = ["$S = 0$"];
    title(title1,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
legend("$Gao\ and\ Wang\ model$","$New\ model$", ...
        'Interpreter', 'latex', 'FontWeight', 'bold','box','off')


    subplot(2,2,3)
    semilogx(Hdd,phase1,color{1},'LineWidth',1)
    hold on
    semilogx(Hdd,phase2,color2{1},'LineWidth',1)
    xlabel('$H_{D}$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$arg(H(i\omega_{D:M2}))$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
        hold on
    plot([1e-8,2e-6],[12.5,12.5],'-b')
 % error
    subplot(2,2,2)
    semilogx(Hdd,100*abs(amplitude1-amplitude2)./amplitude2,color{1},'LineWidth',1)
    hold on
    xlabel('$H_D$','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|\frac{|H_1|-|H_2|}{|H_2|}| \% :M2$','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylim([0,10])
    grid on
    title2 = ["$S = 0$"];
    title(title2,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    
    subplot(2,2,4)
    semilogx(Hdd,100*abs(phase1 - phase2)./abs(phase2),color{1},'LineWidth',1)
    hold on
    xlabel('$H_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}| \% :M2 $','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
    ylim([0,20])
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    hold on
    plot([1e-8,2e-6],[12.5,12.5],'-b')
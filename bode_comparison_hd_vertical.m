%% For tidal behavior of a well in relatively leaky reservoirs
% keep k_overburden is constant
% change the values of h, zw, z
% subroutine needed: solve_tan_eq_n.m
% vertical well
clc;close all;clear
h = [10,60,300];
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

period_components = [period_M2, period_K1, period_Mf, period_Mm, period_Ssa]; %s
nondim_co = kr/(fai*mu*ct*(rw)^2);
period_components_nondim = nondim_co * period_components;
omega_components_non = 2*pi./period_components_nondim;

omega_non  = logspace(-10,-6,400)*sqrt(-1);

S = 1; 

%%%%%%%%%%%%%% semi confined reservoir
color = {'-r','-b','-g','-k'};
color2 = {'--r','--b','--g','--k'};

%% overburden layer
b = 100;
bd = b / rw;
k_overburden = 4e-13;
Hdd = k_overburden / kr / bd ./ hdd;

%% previous model
for jj = 1 : length(Hdd)
    Hd = Hdd(jj);
    hd = hdd(jj);
    C_D = Cdd(jj);

    for i = 1 : length(omega_non)
        omega = omega_non(i) /sqrt(-1) * nondim_co;

        Sd = 1/ (2 * C_D);
        Td = 2 * pi ^2 * kr * h(jj) ./(C(jj) * mu * omega);
        ad = sqrt(2 * pi * sqrt(-1) * Sd/Td);
        betad = sqrt(Hd + ad^2);

        H = (ad/betad)^2 / (1 + ad^2/ (2 * Sd * betad) * besselk(0, betad) / besselk(1, betad) ...
            + pi * sqrt(-1) * S /Td);

        amplitude1(i) = double(abs(H));
        phase1(i) =  double(angle(H))/pi * 180;
    end
%% new model
    lambda_d = Hd * hd;

    x = solve_tan_eq_n(hd*lambda_d,200);
    xx = x/hd;

    for i = 1 : length(omega_non)
        
        s = omega_non(i);

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
        amplitude2(i) = double(abs(H));
        phase2(i) =  double(angle(H))/pi * 180;
    end

    
    fig = figure(1);
    fig.Position = [100 100 800 600];
    subplot(2,1,1)
    semilogx(omega_non/sqrt(-1),amplitude1,color{jj},'LineWidth',1)
    hold on
    semilogx(omega_non/sqrt(-1),amplitude2,color2{jj},'LineWidth',1)
    xlabel('$\omega_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|H(i\omega_D)|$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
  
    grid on

    title1 = ["$S = 1, C_D ="+ string(C_D)+"$"];
    title(title1,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';



    subplot(2,1,2)
    semilogx(omega_non/sqrt(-1),phase1,color{jj},'LineWidth',1)
    hold on
    semilogx(omega_non/sqrt(-1),phase2,color2{jj},'LineWidth',1)
    xlabel('$\omega_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$arg(H(i\omega_D))$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
 % error
    fig = figure(2);
    fig.Position = [100 100 800 600];
    subplot(2,1,1)
    semilogx(omega_non/sqrt(-1),100*abs(amplitude1-amplitude2)./amplitude2,color{jj},'LineWidth',1)
    hold on
    xlabel('$\omega_D$','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|\frac{|H_1|-|H_2|}{|H_2|}| \%$','FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylim([0,20])
    grid on
    title2 = ["$S = 1, C_D ="+ string(C_D)+"$"];
    title(title2,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    
    subplot(2,1,2)
    semilogx(omega_non/sqrt(-1),100*abs(phase1 - phase2)./abs(phase2),color{jj},'LineWidth',1)
    hold on
    xlabel('$\omega_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}| \% $','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
    ylim([0,20])
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
end

tides_name = {'$M_2$','$K_1$','$M_f$','$M_m$','$S_{sa}$'};


figure(1)
text_y = 50;
subplot(2,1,1)
for i = 1 : length(omega_components_non)
    hold on
    plot([omega_components_non(i),omega_components_non(i)],[0,1],'bo-','LineWidth',1)

end
criterion = hdd.*sqrt(Hdd);
legend("$h_D\sqrt{H_D}="+sprintf('%0.3f',criterion(1))+"-Gao\ model$", ...
    "$h_D\sqrt{H_D}="+sprintf('%0.3f',criterion(1))+"-new\ model$", ...
        "$h_D\sqrt{H_D}="+sprintf('%0.3f',criterion(2))+"-Gao\ model$", ...
        "$h_D\sqrt{H_D}="+sprintf('%0.3f',criterion(2))+"-new\ model$", ...
         "$h_D\sqrt{H_D}="+sprintf('%0.3f',criterion(3))+"-Gao\ model$", ...
         "$h_D\sqrt{H_D}="+sprintf('%0.3f',criterion(3))+"-new\ model$", ...
        'Interpreter', 'latex', 'FontWeight', 'bold','box','off')
subplot(2,1,2)
for i = 1 : length(omega_components_non)
    hold on
    plot([omega_components_non(i),omega_components_non(i)],[0,100],'bo-','LineWidth',1)
    text(omega_components_non(i), text_y,tides_name{i},'FontSize',18,'Interpreter','latex','Color','b')
end   

figure(2)
text_y = 10;
subplot(2,1,1)
for i = 1 : length(omega_components_non)
    hold on
    plot([omega_components_non(i),omega_components_non(i)],[0,20],'bo-','LineWidth',1)
    text(omega_components_non(i), text_y,tides_name{i},'FontSize',18,'Interpreter','latex','Color','b')

end
subplot(2,1,2)
for i = 1 : length(omega_components_non)
    hold on
    plot([omega_components_non(i),omega_components_non(i)],[0,20],'bo-','LineWidth',1)
end   
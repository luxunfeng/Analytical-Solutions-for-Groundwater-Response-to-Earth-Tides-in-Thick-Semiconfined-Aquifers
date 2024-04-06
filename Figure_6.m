%% For tidal behavior of a well in relatively leaky reservoirs
% generate Figure 6
% keep k_overburden/k/b is constant 
% subroutine needed: solve_tan_eq_n.m
% vertical well
clc;close all;clear
h = [10:20:600];
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

period_components = [period_M2, period_K1, period_Mf, period_Mm, period_Ssa]; %s
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
fig.Position = [100 100 1600 1200];
for i = 1 : length(omega_components_non)
    for ii = 1 : length(S)
            for jj = 1 : length(hdd)
                Hd = Hdd(jj);
                hd = hdd(jj);
                C_D = Cdd(jj);
            
                    omega = omega_components_non(i) /sqrt(-1) * nondim_co;
            
                    Sd = 1/ (2 * C_D);
                    Td = 2 * pi ^2 * kr * h(jj) ./(C(jj) * mu * omega);
                    ad = sqrt(2 * pi * sqrt(-1) * Sd/Td);
                    betad = sqrt(Hd + ad^2);
            
                    H = (ad/betad)^2 / (1 + ad^2/ (2 * Sd * betad) * besselk(0, betad) / besselk(1, betad) ...
                        + pi * sqrt(-1) * S(ii) /Td);
            
                    amplitude1(ii,jj) = double(abs(H));
                    phase1(ii,jj) =  double(angle(H))/pi * 180;
           
        %% new model
                lambda_d = Hd * hd;
            
                x = solve_tan_eq_n(hd*lambda_d,200);
                xx = x/hd;
                    
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
              
                H = (sum2 + 1)/(1 + C_D * S(ii) * s + C_D * s * sum1);
                amplitude2(ii,jj) = double(abs(H));
                phase2(ii,jj) =  double(angle(H))/pi * 180;
                clear lambda_d x xx s sum1 sum2 H
            end
    end
    tides_name_R = {'$M_2:|\frac{|H_1|-|H_2|}{|H_2|}|$',
        '$K_1:|\frac{|H_1|-|H_2|}{|H_2|}|$',
        '$M_f:|\frac{|H_1|-|H_2|}{|H_2|}|$',
        '$M_m:|\frac{|H_1|-|H_2|}{|H_2|}|$',
        '$S_{sa}:|\frac{|H_1|-|H_2|}{|H_2|}|$'};
    tides_name_td = {'$M_2:|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}|$',
        '$K_1:|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}|$',
        '$M_f:|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}|$',
        '$M_m:|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}|$',
        '$S_{sa}:|\frac{arg(H_1)-arg(H_2)}{arg(H_2)}|$'};
    


    subplot(2,5,i)
    [S_grid, hdd_grid] = meshgrid(S, hdd);
    % here, it is criterion
    % [S_grid, hdd_grid] = meshgrid(S, criterion);
    S_grid = S_grid';
    hdd_grid = hdd_grid';
    error_R = abs(amplitude2-amplitude1)./abs(amplitude2);
    [Cc,hc]=contour(S_grid,hdd_grid,error_R, ...
        'k-','LineWidth',1,'ShowText','on','LabelSpacing',200);
    clabel(Cc, hc, 'FontSize', 16, 'Color', 'k', 'FontWeight', 'normal','interpreter','latex');
    xlabel('$S$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$h_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')

    grid on
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    title1 = tides_name_R{i};
    title(title1,'FontSize',16,'interpreter','latex', 'Color','b','FontWeight', 'bold');
   subplot(2,5,5+i)


   error_td = abs(phase2-phase1)./abs(phase2);
   if i ==1 
       [Cc,hc]=contour(S_grid,hdd_grid,error_td, ...
       [0.02,0.05,0.1,0.2,0.5,1,2,4],'k-', ...
           'LineWidth',1,'ShowText','on','LabelSpacing', 200);
       clabel(Cc, hc, 'FontSize', 16, 'Color', 'k', 'FontWeight', 'normal','interpreter','latex');
   elseif i ==2
        [Cc,hc]=contour(S_grid,hdd_grid,error_td, ...
        [0.01,0.02,0.05,0.1,0.2,0.5,0.7],'k-', ...
           'LineWidth',1,'ShowText','on','LabelSpacing', 200);
       clabel(Cc, hc, 'FontSize', 16, 'Color', 'k', 'FontWeight', 'normal','interpreter','latex');
   else
        [Cc,hc]=contour(S_grid,hdd_grid,error_td, ...
           'k-', ...
           'LineWidth',1,'ShowText','on','LabelSpacing', 200);
       clabel(Cc, hc, 'FontSize', 16, 'Color', 'k', 'FontWeight', 'normal','interpreter','latex');
   end
    
    xlabel('$S$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$h_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
    ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    title1 = tides_name_td{i};
    title(title1,'FontSize',16,'interpreter','latex', 'Color','b','FontWeight', 'bold');
    
end



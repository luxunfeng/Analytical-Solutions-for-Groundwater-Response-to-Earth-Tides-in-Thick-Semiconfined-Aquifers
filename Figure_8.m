%% For tidal behavior of a well in relatively leaky reservoirs
% keep k_overburden/k/b is constant
% subroutine needed: solve_tan_eq_n.m
% Horizontal well
clc;close all;clear
%% parameters
h = [10:20:600];
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

period_components = [period_M2, period_K1, period_Mf, period_Mm, period_Ssa]; %s
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
fig.Position = [100 100 1600 1200];
%% previous model
for i = 1 : length(omega_components_non)
    for ii = 1 : length(S)
        for jj = 1 : length(Hdd)
            Hd = Hdd(jj);
            hd = hdd(jj);
            zd = zdd(jj);
            zwd = zwdd(jj);
        
            
            s = omega_components_non(i);
            sum = 0;
            for j = 1 : 120
                sum = sum + asinh(pi/2/sqrt(s + Hd + (j*pi/hd)^2))*cos(j*pi*zd/hd)*cos(j*pi*zwd/hd);
            end
            A_bar =4/s * 1/ hd/2 *asinh(pi/2/sqrt(s + Hd)) +  8/s * 1 / hd / 2 * sum;
            A = s /(s + Hd);
            H = A/(1 + C_D * S(ii) * s + C_D * s^2 * A_bar);
            amplitude1(ii,jj) = double(abs(H));
            phase1(ii,jj) =  double(angle(H))/pi * 180;
        
        
    %% new modelg 
            lambda_d = Hd * hd;
            x = solve_tan_eq_n(hd*lambda_d,200);
            xx = x/hd;
            s = omega_components_non(i);
            sum1 = 0;
            for j = 1 : 200
                if lambda_d == 0 & j==1
                    sum1 = sum1 + asinh(pi/2/sqrt(s))* 1/2/hd;
                else
                    sum1 = sum1 + asinh(pi/2/sqrt(s  + (xx(j))^2))*cos(xx(j)*zd)*cos(xx(j)*zwd)...
                    * (xx(j)^2 + lambda_d^2)/(hd * (xx(j)^2 + lambda_d^2) + lambda_d);
                end
            end
            sum1 = 4 * sum1;
            
            sum2 = 0;
            for j = 1 : 200
                if lambda_d ==0 
                    sum2 = 0;
                else 
                    sum2 = sum2 - lambda_d*2 / (s + xx(j)^2)*cos(xx(j)*zd)*cos(xx(j)*hd)...
                    * (xx(j)^2 + lambda_d^2)/(hd * (xx(j)^2 + lambda_d^2) + lambda_d);
                end
            end
          
            H = (sum2 + 1)/(1 + C_D * S(ii) * s + C_D * s * sum1);
            amplitude2(ii,jj) = double(abs(H));
            phase2(ii,jj) =  double(angle(H))/pi * 180;
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


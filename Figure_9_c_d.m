%% For tidal behavior of a well in relatively leaky reservoirs
% generate Figure 9 (c) and (d)
% k_overburden changes
% compute the flow rate distribution along the wellbore variation with k'
% subroutine needed: solve_tan_eq_n.m
% the parameters are from wang's paper
% vertical well
clc;close all;clear
h = 48;
rho = 1000;
g = 10;
S = 2.7E-5;
T = 9.6e-6;
fai_ct = S/(h*rho*g);


rw=0.11;
rc = 0.0365;
xw = 0; x = rw;
yw = 0; y = 0;

z = 0:2:h; % the location of measurement gauge
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
% b = 27.70;
% b = 100.00;
 b = 8;
bd = b / rw;
% k_overburden = [0.01:0.04:2]*kr;
k_overburden = [0.005:0.005:0.25]*kr;
Hdd = k_overburden / kr / bd ./ hdd;
kratio = k_overburden/kr;
criteria = hdd.*sqrt(Hdd);
[minvalue_c, minvalue_location_c] = min(abs(criteria-0.2449));
[minvalue_kratio, minvalue_location_kratio] = min(abs(k_overburden / kr-0.1));
max_criteria = max(criteria);
max_kratio = max(kratio);

for jj = 1 : length(Hdd)

    Hd = Hdd(jj);
    hd = hdd(1);
    C_D = Cdd(1);

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

    H = (sum2 + 1)./(1 + C_D * S * s + C_D * s * sum1);
    amplitude2(:,jj) = double(abs(H));
    phase2(:,jj) =  double(angle(H))/pi * 180;


end

fig = figure(1);
fig.Position = [100 100 800 400];
color = colormap("jet");
subplot(1,2,1)
for i = 1:length(Hdd)
    if i == minvalue_location_c
        plot(amplitude2(:,i),zd,'LineWidth',2,'Color','black')
        % plot(amplitude2(:,i),zd,'LineWidth',1,'Color',color(5*i,:))
        txt = text(amplitude2(end,i),zd(end),['$h_D\sqrt{H_D}='+string(criteria(i))+'$'], ...
            "FontSize",18,'Interpreter','latex');
        set(txt, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom');
    elseif i == minvalue_location_kratio
        % plot(amplitude2(:,i),zd,'LineWidth',2,'Color','red')
        plot(amplitude2(:,i),zd,'LineWidth',1,'Color',color(5*i,:))
        % txt = text(amplitude2(end,i),zd(end),['${k^\prime}/{k}='+string(kratio(i))+'$'], ...
        %     "FontSize",18,'Interpreter','latex','Color','red');
        % set(txt, 'HorizontalAlignment', 'center', ...
        %     'VerticalAlignment', 'top');
    else
    plot(amplitude2(:,i),zd,'LineWidth',1,'Color',color(5*i,:))
    hold on
    end
end
xlabel('$|H(i\omega_{D:M2})|$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('$z_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')

grid on

title1 = ["$\frac{h_D}{b_D}="+string(hd/bd)+"$"];
title(title1,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
ax = gca;
set(ax, 'FontSize', 18);
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
cb=colorbar;
% caxis([min(k_overburden/kr),max(k_overburden/kr)])
caxis([min(criteria),max(criteria)])
set(cb, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'LineWidth', 0.5);
ylabel(cb,'$h_D\sqrt{H_D}$','FontSize',18,'Interpreter','latex')
cb.Label.Position = [1.1, 1.07*max_criteria, 0];
cb.Label.Rotation = 0; % set rotation angle



subplot(1,2,2)
for i = 1:length(Hdd)
    plot(phase2(:,i),zd,'LineWidth',1,'Color',color(5*i,:))
    hold on
end
xlabel('$arg(H(i\omega_{D:M2}))$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
% ylabel('$z_D$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')

grid on

title1 = ["$\frac{h_D}{b_D}="+string(hd/bd)+"$"];
title(title1,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
ax = gca;
set(ax, 'FontSize', 18);
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
cb=colorbar;
% caxis([min(k_overburden/kr),max(k_overburden/kr)])
% set(cb, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'LineWidth', 0.5);
% ylabel(cb,'${k^\prime}/{k}$','FontSize',18,'Interpreter','latex')
% cb.Label.Position = [1.1, 1.07*max_kratio, 0];
% cb.Label.Rotation = 0; % set rotation angle
caxis([min(criteria),max(criteria)])
set(cb, 'FontSize', 14, 'TickLabelInterpreter', 'latex', 'LineWidth', 0.5);
ylabel(cb,'$h_D\sqrt{H_D}$','FontSize',18,'Interpreter','latex')
cb.Label.Position = [1.1, 1.07*max_criteria, 0];
cb.Label.Rotation = 0; % set rotation angle
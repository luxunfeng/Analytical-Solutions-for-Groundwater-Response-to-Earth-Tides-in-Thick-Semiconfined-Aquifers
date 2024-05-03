% find another criterion:
% Examine the assumption: uniform flow rate along the wellbore
clear; close all; clc
format long
z = [0.1:0.01:1];

fai= 0.25;
ct = 1e-9;
fai_ct = ct*fai;


rw=0.11;
rc = 0.0365;
xw = 0; x = rw;
yw = 0; y = 0;

mu=0.001;
kr= 1e-13;
kx=kr;ky=kr;kz=kr;
% M2 period
% period_M2 = 12.421*3600; % change frequency here
% period_M2 = 12.421*3600;
% period_K1 = 23.934*3600;
period_Mf = 13.661*24*3600;
% period_Mm = 27.555*24*3600;
% period_Ssa = 0.5*365*24*3600;

period_components = [period_Mf]; %s
nondim_co = kr/(fai_ct*mu*(rw)^2);
period_components_nondim = nondim_co * period_components;
omega_components_non = 2*pi./period_components_nondim*sqrt(-1);

% examine hd effect
fig = figure(1)
fig.Position = [100 100 400 600];
hd = [1:1:100];

for m = 1 : length(hd)
    for i = 1 : length(z)
        x = solve_tan_eq_n(z(i)^2,400);

        sum1 = 0;
        sum2 = 0;
        for j = 1 : 200
            % we still need to consider the effect of hd, this means hd^2*s
            sum1 = sum1 + (x(j)^2 + z(i)^4)/(x(j)^2 + z(i)^4 + z(i)^2) ...
                *z(i)^2/(x(j)^2+(hd(m))^2*omega_components_non) * cos(x(j))^2;

            sum2 = sum2 + (x(j)^2 + z(i)^4)/(x(j)^2 + z(i)^4 + z(i)^2) ...
                *z(i)^2/(x(j)^2+(hd(m))^2*omega_components_non) * cos(x(j));
        end
        sum1 = sum1 * 2;
        sum2 = sum2 * 2;
        ratio(i) = abs((1-sum1)/(1-sum2));

    end
    if m == length(hd)
        plot(z,ratio,'-g','LineWidth',2)
    end

    plot(z,ratio)
    hold on


end
xlabel('$h_D\sqrt{H_D}$','Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('$\frac{|H:{z_D=h_D}|}{|H:{z_D=0}|},\frac{|p_{wD}:z_D=h_D|}{|p_{wD}:z_D=0|}$','Interpreter', 'latex', 'FontWeight', 'bold')
  
ax = gca;
    set(ax, 'FontSize', 18); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
% assume h^2/k > 1e13
omega = 2 * pi / period_components;
h2k = 1e+13;
for i = 1 : length(z)
    x = solve_tan_eq_n(z(i)^2,400);

    sum1 = 0;
    sum2 = 0;
    for j = 1 : 200
        % we still need to consider the effect of hd, this means hd^2*s
        sum1 = sum1 + (x(j)^2 + z(i)^4)/(x(j)^2 + z(i)^4 + z(i)^2) ...
            *z(i)^2/(x(j)^2+omega*h2k*fai_ct*mu* ...
            sqrt(-1)) * cos(x(j))^2;

        sum2 = sum2 + (x(j)^2 + z(i)^4)/(x(j)^2 + z(i)^4 + z(i)^2) ...
            *z(i)^2/(x(j)^2+omega*h2k*fai_ct*mu* ...
            sqrt(-1)) * cos(x(j));
    end
    sum1 = sum1 * 2;
    sum2 = sum2 * 2;
    ratio(i) = abs((1-sum1))/abs((1-sum2));

end
hold on
plot(z,ratio,'-b','LineWidth',2)
hold on
plot([0.1,0.475],[0.9,0.9],'-k','LineWidth',1.5)
title('$M_f$','Interpreter','latex','FontSize',16,'FontWeight','bold','FontAngle','italic')
xlim([0.1,1])
ylim([0.6,1.1])
grid on
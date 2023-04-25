% find another criterion:
% a = hd*sqrt(HD)
% hd and omega also have effects.
% here it is hd's effect
clear; close all; clc
format long
z = [0.1:0.01:1];
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

mu=0.001;
kr=T/h*mu/rho/g;
kx=kr;ky=kr;kz=kr;
% M2 period
period_M2 = 12.421*3600; % change frequency here


period_components = [period_M2]; %s
nondim_co = kr/(fai_ct*mu*(rw)^2);
period_components_nondim = nondim_co * period_components;
omega_components_non = 2*pi./period_components_nondim*sqrt(-1);

% examine hd effect
figure

hd = [1:1:100];
for m = 1 : length(hd)
    for i = 1 : length(z)
        x = solve_tan_eq_n(z(i)^2,400);

        sum1 = 0;
        sum2 = 0;
        for j = 1 : 200
            % we still need to consider the effect of hd, this means hd^2*s
            sum1 = sum1 + (x(j)^2 + z(i)^4)/(x(j)^2 + z(i)^4 + z(i)^2) ...
                *z(i)^2/(x(j)^2+(hd(m))^2*omega_components_non* ...
                sqrt(-1)) * cos(x(j))^2;
            %         sum2 = sum2 + a(i)^2/x(j)^2 * cos(x(j));
            sum2 = sum2 + (x(j)^2 + z(i)^4)/(x(j)^2 + z(i)^4 + z(i)^2) ...
                *z(i)^2/(x(j)^2+(hd(m))^2*omega_components_non* ...
                sqrt(-1)) * cos(x(j));
        end
        sum1 = sum1 * 2;
        sum2 = sum2 * 2;
        ratio(i) = abs((1-sum1)/(1-sum2));

    end
plot(z,ratio)
hold on

end
xlabel('criterion')
ylabel('pressure ratio')
ylim([0,1])

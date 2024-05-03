clear all; clc; close all;
period_m2 = 12.421; % hours
frequency_m2 = 2*pi/period_m2/3600; % rad/s

T_series = [1e-4,1e-6,1e-8];
% T_series = [1e-8];
S_series = [1e-4,1e-6,1e-8];
% color_series = {'r','b','g','k','c','m'};
% marker_series = {'o','x','*','+'};

marker_size = 6;
line_width = 1.5;
marker_interval = 5;
figure(1)

rc = 0.1;
rw = 0.1;
K_b  = 10.^[-10, -9, -8, -7, -5, -4];
%% HERBERT WANG
% for j = 1 : length(T_series)
for j = 1
    % for k = 1 : length(S_series)
    for k = 2
        h = 48;
        rho = 1000;
        g = 10;
        fai_ct = S_series(k)/ (h*rho*g); % this is from Xuhua's thesis or infered by comparison between governing equations.
        xw = 0; x = rw;
        yw = 0; y = 0;
        mu=0.001;
        kr = T_series(j) /h * mu / rho / g;
        kx = kr; ky = kr; kz = kr;
        etax=kx/(fai_ct*mu);etaz=kz/(fai_ct*mu);etay=ky/(fai_ct*mu);

        lambda_wang = h * sqrt(i*frequency_m2/etax);
        z = 0:0.001:h;
        H = 1 + tanh(lambda_wang)*sinh(lambda_wang*z/h) - cosh(lambda_wang*z/h);
        figure(1)
        subplot(1,2,1);
        hold on
        % yyaxis left
        plot(abs(H),(h-z)/h, 'k-','LineWidth',line_width);
        % set(gca, 'YDir', 'reverse')
        subplot(1,2,2);
        hold on
        % yyaxis left
        plot(angle(H)/pi*180,(h-z)/h, 'k-','LineWidth',line_width);
        % set(gca, 'YDir', 'reverse')
    end
end

% %% Wang 2018
% % for j = 1 : length(T_series)
% for j = 1
%     % for k = 1 : length(S_series)
%     for k = 3
%         h = 48;
%         rho = 1000;
%         g = 10;
%         fai_ct = S_series(k)/ (h*rho*g); % this is from Xuhua's thesis or infered by comparison between governing equations.
%         xw = 0; x = rw;
%         yw = 0; y = 0;
%         mu=0.001;
%         kr = T_series(j) /h * mu / rho / g;
%         kx = kr; ky = kr; kz = kr;
%         etax=kx/(fai_ct*mu);etaz=kz/(fai_ct*mu);etay=ky/(fai_ct*mu);
%         period_M2 = 12.421*3600;
%         period_components = [period_M2]; %s
%         nondim_co = kr/(fai_ct*mu*(rw)^2);
%         period_components_nondim = nondim_co * period_components;
%         omega_components_non = 2*pi./period_components_nondim*sqrt(-1);
% 
%        H_D = K_b/T_series(j) * rw^2;
%        H = omega_components_non ./ (H_D + omega_components_non);
%        for m = 1 : length(K_b)
%         figure(1)
%         subplot(1,2,1);
%         hold on
%         % yyaxis left
%         plot(abs(H(m))*ones(length(z),1),z/h, '-','Color',color_series{m},'LineWidth',line_width);
%         % set(gca, 'YDir', 'reverse')
%         subplot(1,2,2);
%         hold on
%         % yyaxis left
%         plot(angle(H(m))/pi*180 * ones(length(z),1),z/h, '-','Color',color_series{m},'LineWidth',line_width);
%         % set(gca, 'YDir', 'reverse')
%        end
%     end
% end


%% LSH model
h = 48;
rho = 1000;
g = 10;
% T_series
% for j = 1 : length(T_series)
for j = 1
    % for k = 1 : length(S_series)
    for k = 2
        fai_ct = S_series(k)/ (h*rho*g); % this is from Xuhua's thesis or infered by comparison between governing equations.
        xw = 0; x = rw;
        yw = 0; y = 0;
        z_series = 0:2:h; % the location of measurement gauge

        mu=0.001;
        kr = T_series(j) /h * mu / rho / g;
        kx = kr; ky = kr; kz = kr;

        etax=kx/(fai_ct*mu);etaz=kz/(fai_ct*mu);etay=ky/(fai_ct*mu);

        C = pi * rc^2/ rho / g; % m3/Pa, from defination of wellbore storage coeffi
        C=0; % here, C should be set to zero to filter the effect of wellbore flow.

        CD = C./(2*pi*h*fai_ct*rw^2); %non-dimensional wellbore storage

        rd = rw / rw;
        hd=h/(rw);


        period_M2 = 12.421*3600;
        period_components = [period_M2]; %s
        nondim_co = kr/(fai_ct*mu*(rw)^2);
        period_components_nondim = nondim_co * period_components;
        omega_components_non = 2*pi./period_components_nondim*sqrt(-1);
        Skin = 0;
        lambda_d = K_b * rw * mu / rho / g / kr;
        B_LSH = sqrt(hd*lambda_d);

        for m = 1 : length(lambda_d)
            for kk = 1 : length(z_series)
                z = z_series(kk); % it is noted that the z direction is opposite in our study (upwards) compared with Herbert Wang
                zd = z/(rw);
                if hd*lambda_d(m) > 1e+6
                    tol = 1e-10;
                    x = solve_tan_eq_n(hd*lambda_d(m),40000,tol);
                else
                    x = solve_tan_eq_n(hd*lambda_d(m),40000);
                end

                xx = x/hd;
                s = omega_components_non(1);
                sum1 = 0;
                for jj = 1 : length(x)
                    sum1 = sum1 + 1/xx(jj) * sin(xx(jj) * hd) * cos(xx(jj)*zd)...
                        * besselk(0, sqrt(s + xx(jj)^2))...
                        * (xx(jj)^2 + lambda_d(m)^2)/(hd * (xx(jj)^2 + lambda_d(m)^2) + lambda_d(m));
                end
                sum1 = 2 * sum1;
                sum2 = 0;
                for jj = 1 : length(x)
                    sum2 = sum2 - lambda_d(m)*2/ (s + xx(jj)^2)*cos(xx(jj)*zd)*cos(xx(jj)*hd)...
                        * (xx(jj)^2 + lambda_d(m)^2)/(hd * (xx(jj)^2 + lambda_d(m)^2) + lambda_d(m));
                end

                H = (sum2 + 1)/(1 + CD * Skin * s + CD * s * sum1);
                amplitude2(kk) = double(abs(H));
                phase2(kk) =  double(angle(H))/pi * 180;
            end
            m
            figure(1)
            subplot(1,2,1);
            % yyaxis right
            % plot(amplitude2,(z_series)/h, 'o','LineWidth',line_width,'MarkerSize',m,'Color',color_series{m})
            plot(amplitude2,(z_series)/h, 'o','LineWidth',line_width,'MarkerSize',m)
            hold on
    xlabel('$Amplitude\ ratio$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$z/h $','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
        ax = gca;
    set(ax, 'FontSize', 20); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
            % set(gca, 'YDir', 'reverse')

            subplot(1,2,2)
            % yyaxis right
            % plot(phase2,(z_series)/h , 'o','LineWidth',line_width,'MarkerSize',m,'Color',color_series{m})
            plot(phase2,(z_series)/h , 'o','LineWidth',line_width,'MarkerSize',m)

            hold on
    xlabel('$Phase\ shift$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    ylabel('$z/h $','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold')
    grid on
        ax = gca;
    set(ax, 'FontSize', 20); 
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
            % set(gca, 'YDir', 'reverse')
        end
    end
end

% lg =legend('$Detournay \& Cheng, 1993$','$New model\ \frac{K^\prime}{b^\prime}=10^{-10}$', '$New model\ \frac{K^\prime}{b^\prime}=10^{-9}$','$New model\ \frac{K^\prime}{b^\prime}=10^{-8}$','$New model\ \frac{K^\prime}{b^\prime}=10^{-7}$','$New model\ \frac{K^\prime}{b^\prime}=10^{-5}$','$New model\ \frac{K^\prime}{b^\prime}=10^{-4}$',...
% 'FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold','box','off')
% lg.Location = 'eastoutside'
% function for plot

function semilog_with_markers(x, y, lineStyle, markerStyle, markerInterval,line_width,marker_size)
semilogx(x, y, lineStyle,'LineWidth',line_width);      % Plot the basic line
hold on;
markerIndices = 1:markerInterval:length(x); % Calculate indices for markers
semilogx(x(markerIndices), y(markerIndices), markerStyle,'MarkerSize',marker_size); % Plot markers
hold off;
end
function loglog_with_markers(x, y, lineStyle, markerStyle, markerInterval, line_width, marker_size)
loglog(x, y, lineStyle,'LineWidth',line_width);      % Plot the basic line
hold on;
markerIndices = 1:markerInterval:length(x); % Calculate indices for markers
loglog(x(markerIndices), y(markerIndices), markerStyle,'MarkerSize',marker_size); % Plot markers
hold off;
end

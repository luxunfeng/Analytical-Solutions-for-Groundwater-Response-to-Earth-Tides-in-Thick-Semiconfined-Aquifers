clear all; clc; close all;
period_m2 = 12.421; % hours
frequency_m2 = 2*pi/period_m2/3600; % rad/s


T_series = [1e-4,1e-6,1e-8];
% T_series = [1e-8];
S_series = [1e-4,1e-6,1e-8];
color_series = {'r','b','g','k'};
% marker_series = {'o','x','*','+'};

marker_size = 6;
line_width = 1.5;
marker_interval = 5;
figure(1)
for j = 1 : length(T_series)
    for k = 1 : length(S_series)
        T = T_series(j);
        S = S_series(k);
        rc = 0.1;
        rw = 0.1;
        K_b  = 10.^[-15:0.1:-3];

        beta = (1/T * K_b + i * frequency_m2 * S / T).^0.5;

        epsilon = 1 + (rc/rw)^2 * i * frequency_m2 * rw / 2 / T ./beta .* besselk(0, beta * rw)./...
            besselk(1, beta*rw);


        H = i * frequency_m2 * S./ (i * frequency_m2 * S + K_b) ./ epsilon;

        amplitude = abs(H);
        phase = angle(H)/pi * 180;

        subplot(2,length(T_series),j)
        % loglog_with_markers(K_b, amplitude,['-',color_series{k}],[color_series{k},marker_series{k}], ...
        % marker_interval,line_width,marker_size)
        semilogx(K_b, amplitude, ['-',color_series{k}],'LineWidth',line_width)
        hold on
        xlim([min(K_b),max(K_b)])
        ylim([0,1])

        subplot(2,length(T_series),j + 3)
        % semilog_with_markers(K_b, phase,['-',color_series{k}],[color_series{k},marker_series{k}], ...
        % marker_interval,line_width,marker_size)
        semilogx(K_b, phase, ['-',color_series{k}],'LineWidth',line_width)
        hold on
        xlim([min(K_b),max(K_b)])
        ylim([-90,90])


    end
end



%% LSH model

h = 48;
rho = 1000;
g = 10;
% T_series
for j = 1 : length(T_series)
% for j = 1
    for k = 1 : length(S_series)
    % for k = 3
        fai_ct = S_series(k)/ (h*rho*g); % this is from Xuhua's thesis or infered by comparison between governing equations.
        xw = 0; x = rw;
        yw = 0; y = 0;
        z = 0; % the location of measurement gauge
        mu=0.001;
        kr = T_series(j) /h * mu / rho / g;
        kx = kr; ky = kr; kz = kr;

        etax=kx/(fai_ct*mu);etaz=kz/(fai_ct*mu);etay=ky/(fai_ct*mu);

        C = pi * rc^2/ rho / g; % m3/Pa, from defination of wellbore storage coeffi
        % C=0;

        % h = 0.4*h;
        CD = C./(2*pi*h*fai_ct*rw^2); %non-dimensional wellbore storage

        rd = rw / rw;
        hd=h/(rw);
        zd = z/(rw);

        period_M2 = 12.421*3600;
        period_components = [period_M2]; %s
        nondim_co = kr/(fai_ct*mu*(rw)^2);
        period_components_nondim = nondim_co * period_components;
        omega_components_non = 2*pi./period_components_nondim*sqrt(-1);
        Skin = 0;
        lambda_d = K_b * rw * mu / rho / g / kr;
        B_LSH = sqrt(hd*lambda_d);

        for m = 1 : length(lambda_d)
            if hd*lambda_d(m) > 1e+6
                tol = 1e-10;
                x = solve_tan_eq_n(hd*lambda_d(m),50000,tol);
            else
                x = solve_tan_eq_n(hd*lambda_d(m),50000);
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
            % sum1 = 0;
            sum2 = 0;
            for jj = 1 : length(x)
                sum2 = sum2 - lambda_d(m)*2/ (s + xx(jj)^2)*cos(xx(jj)*zd)*cos(xx(jj)*hd)...
                    * (xx(jj)^2 + lambda_d(m)^2)/(hd * (xx(jj)^2 + lambda_d(m)^2) + lambda_d(m));
            end
            cos(xx(1)*zd);
            sum2;
            Sum2(m) = sum2;
            Sum1(m) = sum1;
            % figure(2)
            %  ss = linspace(100,10,length(xx));
            % %  % c = linspace(1,100,length(xx));
            % %  % c(1) = 10;
            %  ss(1) = 150;
            % % 
            % % % plot(cos(xx*zd).*cos(xx*hd),'o'); hold on
            % % scatter(real(- lambda_d(m)*2./ (s + xx.^2).*cos(xx*zd).*cos(xx*hd)...
            % %         .* (xx.^2 + lambda_d(m)^2)./(hd * (xx.^2 + lambda_d(m)^2) + lambda_d(m))), ...
            % %         imag(- lambda_d(m)*2./ (s + xx.^2).*cos(xx*zd).*cos(xx*hd)...
            % %         .* (xx.^2 + lambda_d(m)^2)./(hd * (xx.^2 + lambda_d(m)^2) + lambda_d(m))),ss); hold on
            % % 
            % % scatter(real(sum2),imag(sum2),150-m,'s','filled'); hold on
            % subplot(1,3,1)
            %            semilogx(real(cos(xx*zd).*sin(xx*hd)...
            %         .* (xx.^2 + lambda_d(m)^2)./(xx.*hd.* (xx.^2 + lambda_d(m)^2) + lambda_d(m))),'-bo'); hold on
            %            subplot(1,3,2)
            %         semilogx( real(xx.^2./(s+xx.^2)),'-rs'); hold on
            %                                subplot(1,3,3)
            %         semilogx( imag(xx.^2./(s+xx.^2)),'-rs'); hold on
            % 
            %  real_ = sum(cos(xx*zd).*sin(xx*hd)...
            %         .* (xx.^2 + lambda_d(m)^2)./(xx.*hd.* (xx.^2 + lambda_d(m)^2) + lambda_d(m)) .* xx.^2./(s+xx.^2))
            %  approximate = sum(cos(xx*zd).*sin(xx*hd)...
            %         .* (xx.^2 + lambda_d(m)^2)./(xx.*hd.* (xx.^2 + lambda_d(m)^2) + lambda_d(m)) .* xx(1).^2./(s+xx(1).^2))
            % 
            %  real(real_) <= real(approximate)
            %  abs(imag(real_)) >=abs(imag(approximate)) 
             % sum2 = 0;

            H = (sum2 + 1)/(1 + CD * Skin * s + CD * s * sum1);
            amplitude2(m) = double(abs(H));
            phase2(m) =  double(angle(H))/pi * 180;

        end
        j
        k
        figure(1)
        subplot(2,length(T_series),j);
        % loglog_with_markers(K_b, amplitude2,['--',color_series{k}],[color_series{k},marker_series{k}], ...
        %     marker_interval,line_width,marker_size)
        % semilogx(ax1, B_LSH, amplitude2, ['--',color_series{k}],'LineWidth',line_width)
        % ax1.XAxisLocation = 'top';
     
        semilogx(K_b, amplitude2, ['--',color_series{k}],'LineWidth',line_width)
        hold on
        % xlim([min(B_LSH),max(B_LSH)])
        xlim([min(K_b),max(K_b)])
        ylim([0,1])
        

        subplot(2,length(T_series),j+3)
        % semilog_with_markers(K_b, phase2,['--',color_series{k}],[color_series{k},marker_series{k}], ...
        %     marker_interval,line_width,marker_size)
        semilogx(K_b, phase2, ['--',color_series{k}],'LineWidth',line_width)
        hold on
        xlim([min(K_b),max(K_b)])
                ylim([-90,90])

        

    end
end

lg =legend('$Storativity=10^{-4}-Wang\ et\ al., 2018$','$Storativity=10^{-6}-Wang\ et\ al., 2018$', ...
    '$Storativity=10^{-8}-Wang\ et\ al., 2018$','$Storativity=10^{-4}-new\ model$', ...
    '$Storativity=10^{-6}-new\ model$','$Storativity=10^{-8}-new\ model$', ...
    'FontSize',16,'Interpreter', 'latex', 'FontWeight', 'bold','box','off')
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

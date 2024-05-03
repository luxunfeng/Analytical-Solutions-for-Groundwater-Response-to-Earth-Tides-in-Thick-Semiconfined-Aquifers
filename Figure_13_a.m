clc;close all; clear
h = 48;
b = 277;
k_ = linspace(3e-15, 3e-14,1000);
% k_ = 3e-14;
k = linspace(2e-14, 3e-12,1000);
for i = 1 : length(k)
    for j = 1 : length(k_)
        results(i,j) = sqrt(h*k_(j)/b/k(i));
    end
end

[K_, K] = meshgrid(k_, k);

[M,c] = contour(K_,K, results,'ShowText','on');
% plot(k,results)

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xticks([3e-15,1e-14,2e-14,3e-14])
yticks([2e-14,1e-13,1e-12, 3e-12])

%修饰
grid on
ax = gca;
set(ax, 'FontSize', 14);
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

c.LineWidth = 1;
c.LineStyle = '-';
c.EdgeColor = 'k';
clabel(M,c, 'FontSize', 12, 'FontWeight', 'bold','interpreter','latex');
title2 = ["$h=48m, b=277m$"];
title(title2,'FontSize',14,'interpreter','latex', 'FontWeight', 'bold');
    xlabel('$Permeability\ of\ overlaying\ layer[m^{2}]$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold','FontAngle','normal')
    ylabel('$Permeability\ of\ target\ reservoir[m^{2}]$','FontSize',14,'Interpreter', 'latex', 'FontWeight', 'bold','FontAngle','normal')
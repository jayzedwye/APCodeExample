function [hd] = plot_fun(coeff_total, std_coeff, horizon, color)

y1 = coeff_total(1,1:horizon);
y2 = coeff_total(1,1:horizon)+2*std_coeff(1:horizon);
y3 = coeff_total(1,1:horizon)-2*std_coeff(1:horizon);

x  = linspace(1,horizon,horizon);

% area does not work with NaN

% plot the area
y = [y3' (y2'-y3')];
ha = area(x, y);hold on;
set(ha(1), 'FaceColor', 'none') % this makes the bottom area invisible
set(ha, 'LineStyle', 'none')

% plot the line edges
hb = plot(x, y3, 'LineWidth', 1);hold on;
hc = plot(x, y2, 'LineWidth', 1);hold on;

% plot the point estimates
hd = plot(x, y1, 'LineWidth', 3);hold on;
%he = plot(x, t1, 'LineWidth', 3,'LineStyle','--');hold on;
% set the line and area colors
set(ha(2), 'FaceColor',[0.7 0.7 0.7]);
set(hb, 'Color', [0.9 0.9 0.9]);
set(hc, 'Color', [0.9 0.9 0.9]);
set(hd, 'Color', 'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1 = coeff_total(1,1:horizon);
y2 = coeff_total(1,1:horizon)+1*std_coeff(1:horizon);
y3 = coeff_total(1,1:horizon)-1*std_coeff(1:horizon);



x  = linspace(1,horizon,horizon);
% area does not work with NaN

% plot the area
y = [y3' (y2'-y3')];
ha = area(x, y);hold on;
set(ha(1), 'FaceColor', 'none') % this makes the bottom area invisible
set(ha, 'LineStyle', 'none')

% plot the line edges
hb = plot(x, y3, 'LineWidth', 1);hold on;
hc = plot(x, y2, 'LineWidth', 1);hold on;

% plot the point estimates
hd = plot(x, y1, 'LineWidth', 3);hold on;
%he = plot(x, t1, 'LineWidth', 3,'LineStyle','--');hold on;
% set the line and area colors
set(ha(2), 'FaceColor',[0.5 0.5 0.5]);
set(hb, 'Color', [0.9 0.9 0.9]);
set(hc, 'Color', [0.9 0.9 0.9]);
set(hd, 'Color', color);


axis([1 horizon -1 1.5])

xlabel('Horizon in years')
ylabel('% of Variance')
set(gca, 'Layer', 'top','FontSize',9)
set(gca, 'FontName', 'Times New Roman')

grid

end
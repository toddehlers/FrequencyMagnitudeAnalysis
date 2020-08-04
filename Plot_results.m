% Plot routines, separated by calculation

clear all;
load('results/5years_all');
addpath('routines');


rockfalldata = dlmread('data/rockfall_data_all.txt');
rockfalldata = sortrows(rockfalldata);
time = 5.2;

V = rockfalldata(:,1);

nevent = numel(V);
numbers = [nevent:-1:1]'/time;

darkblue = [31,120,180]/255;
lightblue = [166,206,227]/255;

darkgreen = [51,160,44]/255;
lightgreen = [178,223,138]/255;

report = @(A) [min(A); mean(A); max(A)];

%% Plot adjusted R^2-value of fit as fuction of removing first events

min_event_ro = 0;
max_event_ro = round(nevent*0.85);

%% Plot adjusted R^2-value of fit as fuction of removing first events

figure('Name', 'Finding Rollover',...
    'Units', 'centimeters',...
    'Position', [1,1,18,14])



spmain = subplot(2,1,1);
V_maxplot = 1.5;

yyaxis left
plot(V(1:max_event_ro),gof_arr(:,1),'LineWidth',2.5)
line([V_ro1 V_ro1], [0.8,1.1],'LineStyle','--','Color',0.3*[1 1 1],'LineWidth',2.0)
line([V_ro2 V_ro2], [0.8,1.1],'LineStyle','--','Color',0.3*[1 1 1],'LineWidth',2.0)
text(V_ro1+0.05, 0.993, '1st peak in R^{2} value', 'Rotation', 90, 'HorizontalAlignment', 'right', ...
    'FontSize', 15, 'Color', 0.3*[1 1 1])
text(V_ro2+0.05, 0.993, '2nd peak in R^{2} value', 'Rotation', 90, 'HorizontalAlignment', 'right', ...
    'FontSize', 15, 'Color', 0.3*[1 1 1])


ylabel('adjusted R^2 of nonlinear fit')


yyaxis right
ksline = plot(V(1:max_event_ro),gof_arr(:,2),'LineWidth',2.5);
ylabel('KS statistic')

%calculate values for annotations and guidelines
[~, i_ro1] = min(abs(V-V_ro1));
[~, i_ro2] = min(abs(V-V_ro2));
[~, i_maxplot] = min(abs(V-V_maxplot));
KS_mean1 = mean(gof_arr(i_ro1:i_r2min1,2));
KS_mean2 = mean(gof_arr(i_ro2:i_maxplot,2));

line([0.95*V_ro1 1.05*V(i_r2min1)], [KS_mean1, KS_mean1],'LineStyle',':',...
    'Color',ksline.Color,'LineWidth',2.0)
line([0.95*V_ro2 1.05*V_maxplot], [KS_mean2, KS_mean2],'LineStyle',':','Color',ksline.Color,'LineWidth',2.0)

axmain=gca;

% Plot second axes, percentage of data discarded
spbot = subplot(2,1,2);

datadiscard = 100*[0:numel(V)-1]./numel(V);
plot(V, datadiscard,'LineWidth',2.5, 'Color',0.5*[1 1 1])
xlabel( 'Rollover Volume [m^3]')
botl = legend('% of data discarded');
set(botl, 'Box', 'off');
set(botl, 'FontSize', 12)
set(botl, 'Location', 'southeast')
axbot = gca;

%set properties of axes, labels and ticks...

for ax = [axmain, axbot]
    set(ax,'XTick', [0:0.25:1.5])
    ax.XLim = [0,1.5];
    ax.LineWidth = 2.0;
    ax.FontSize = 16;
end


axmain.YAxis(1).Limits = [0.91,1.005];
axmain.YAxis(1).TickValues = [0.92:0.02:1];
xmain.YAxis(1).TickLabels = {'0.92', '0.94', '0.96', '0.98', '1.00'};
axmain.YAxis(2).Limits = [0.03,0.15];
set(spmain, 'XAxisLocation', 'top');



axbot.YAxis.Limits = [0, 100];


linkaxes([spmain, spbot], 'x')
set(spmain, 'position', [0.125 0.30 0.75 0.62])
set(spbot, 'position', [0.125 0.15 0.75 0.15])

set(axbot, 'YTick', 0:25:100)
set(axbot, 'YTickLabel', {'0%', '', '50%', '', '100%'})

annotation('textarrow', [0.74 0.66], [0.533 0.352],...
    'String',{'minimum','KS-value'},...
    'FontSize',15, 'Color', 0.3*[1 1 1]); 


%print('Figures/Rollover_Limits', '-dpdf', '-r300')

%% create Range of possible rollovers, log spaced

nstep_ro = 25;
ro_range = exp(linspace(log(min(V_ro1,V_ro2)),log(max(V_ro1,V_ro2)),nstep_ro));

v_eval = 1000; %volume of event for return time calculation
t_eval = 100; %time interval for which the maximum recurring volume is calculate (e.g. 100-year-event)
vmax_er = 1e4; %upper volume limit for erosion integration

%create upper and lower limit for erosion and return time
tret_max = zeros(nstep_ro,3);
tret_min = zeros(nstep_ro,3);
erosion_max = zeros(nstep_ro,3);
erosion_min = zeros(nstep_ro,3);
for i = 1:nstep_ro
    for j= 1:3
        fitres.a = a_res(i,j);
        
        %max values
        fitres.b = -b_res(i,j)+berr_res(i,j)+1;
        [tmax, ~, emax] = eval_fit(fitres, v_eval, t_eval, V_ro1, vmax_er);
        tret_max(i,j) = tmax;
        erosion_max(i,j) = emax;
        
        %min values
        fitres.b = -b_res(i,j)-berr_res(i,j)+1;
        [tmin, ~, emin] = eval_fit(fitres, v_eval, t_eval, V_ro1, vmax_er);
        tret_min(i,j) = tmin;
        erosion_min(i,j) = emin;
    end
end
        
%% Rollover dependency figure

figure('Name', 'Rollover Effect - Shaded regions',...
    'Units', 'centimeters',...
    'Position', [1,1,18,31.78])

ax1 = subplot(3,1,1);
hold on
[xpatch, ypatch] = create_patch(ro_range,b_res(:,3)+berr_res(:,3)-1,b_res(:,3)-berr_res(:,3)-1);
patch('XData', xpatch, 'YData', ypatch,...
    'FaceColor',lightblue,'FaceAlpha', 0.5,...
    'LineStyle','none');
plot(ro_range,b_res(:,3)+berr_res(:,3)-1,'Color',darkblue, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,b_res(:,3)-berr_res(:,3)-1,'Color',darkblue, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,b_res(:,3)-1,'Color',darkblue, 'LineStyle','-','LineWidth',2.5);

[xpatch, ypatch] = create_patch(ro_range,b_res(:,2)+berr_res(:,2)-1,b_res(:,2)-berr_res(:,2)-1);
patch('XData', xpatch, 'YData', ypatch,...
    'FaceColor',lightgreen,'FaceAlpha', 0.5,...
    'EdgeColor',darkgreen,'LineWidth',1.5,'LineStyle','none');
plot(ro_range,b_res(:,2)+berr_res(:,2)-1,'Color',darkgreen, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,b_res(:,2)-berr_res(:,2)-1,'Color',darkgreen, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,b_res(:,2)-1,'Color',darkgreen, 'LineStyle','-','LineWidth',2.5);

%set specific axes values
ylabel('scaling exponent b','FontWeight','bold','FontSize',18);
ylim([ 0.5,0.85])
box on
set(ax1,'XAxisLocation', 'top')
xlabel('rollover volume [m^{3}]','FontWeight','bold','FontSize',18);
text(0.28, 0.82, 'a', 'FontSize',20, 'FontWeight', 'bold')


ax2 = subplot(3,1,2);
hold on
%plot range of return times as patch
[xpatch, ypatch] = create_patch(ro_range,tret_max(:,3),tret_min(:,3));
patch('XData', xpatch, 'YData', ypatch,...
    'FaceColor',lightblue,'FaceAlpha', 0.5,...
    'LineStyle','none');
%plot edges and mean line (bold)
plot(ro_range,tret_min(:,3),'Color',darkblue, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,tret_max(:,3),'Color',darkblue, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,tret_res(:,3),'Color',darkblue, 'LineStyle','-','LineWidth',2.5);

[xpatch, ypatch] = create_patch(ro_range,tret_max(:,2),tret_min(:,2));
patch('XData', xpatch, 'YData', ypatch,...
    'FaceColor',lightgreen,'FaceAlpha', 0.5,...
    'LineStyle','none');
%plot edges and mean line (bold)
plot(ro_range,tret_min(:,2),'Color',darkgreen, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,tret_max(:,2),'Color',darkgreen, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,tret_res(:,2),'Color',darkgreen, 'LineStyle','-','LineWidth',2.5);

%set specific axes values
ylabel('1000 m^3 return time [yr]','FontWeight','bold','FontSize',18);
ylim([ 0, 12])
set(ax2,'XTickLabel', [])
text(0.28, 11.0, 'b', 'FontSize',20, 'FontWeight', 'bold')

ax3 = subplot(3,1,3);
hold on
%plot range of erosion rates as patch
[xpatch, ypatch] = create_patch(ro_range,erosion_max(:,3),erosion_min(:,3));
patch('XData', xpatch, 'YData', ypatch*1e3,...
    'FaceColor',lightblue,'FaceAlpha', 0.5,...
    'LineStyle','none');
%plot edges and mean line (bold)
plot(ro_range,erosion_min(:,3)*1e3,'Color',darkblue, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,erosion_max(:,3)*1e3,'Color',darkblue, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,erosion_res(:,3)*1e3,'Color',darkblue, 'LineStyle','-','LineWidth',2.5);

[xpatch, ypatch] = create_patch(ro_range,erosion_max(:,2),erosion_min(:,2));
patch('XData', xpatch, 'YData', ypatch*1e3,...
    'FaceColor',lightgreen,'FaceAlpha', 0.5,...
    'LineStyle','none');
%plot edges and mean line (bold)
plot(ro_range,erosion_min(:,2)*1e3,'Color',darkgreen, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,erosion_max(:,2)*1e3,'Color',darkgreen, 'LineStyle','-','LineWidth',1.5);
plot(ro_range,erosion_res(:,2)*1e3,'Color',darkgreen, 'LineStyle','-','LineWidth',2.5);

%set specific axes values
xlabel('rollover volume [m^{3}]','FontWeight','bold','FontSize',18);
ylabel('long-term erosion [mm/yr]','FontWeight','bold','FontSize',18);
ylim([0.05,0.3])
text(0.28, 0.2786, 'c', 'FontSize',20, 'FontWeight', 'bold')

%set general axes values
linkaxes([ax1, ax2, ax3], 'x')
xlim([0.25,0.95])
for ax = [ax1, ax2, ax3]
    ax.LineWidth = 2.0;
    ax.FontSize = 16;
    ax.Box = 'on';
    ax.XTick = [0.3:0.15:1.1];
    ax.YAxis.TickLabelFormat = '%,.2f';
end

ax1.YTick = [0.5:0.05:0.85];
ax2.YAxis.TickLabelFormat = '%,.1f';

ax1.Position = [0.1323    0.6624    0.8426    0.2728];
ax2.Position = [0.1323    0.3636    0.8426    0.2728];
ax3.Position = [0.1323    0.0648    0.8426    0.2728];

%print('Figures/Rollover_Effect', '-dpdf', '-r300')

%% Additional data for summary figure

%Create extended error envelope

X = logspace(-2,4,100);
Y_errtop = zeros(100,1);
Y_errbottom = zeros(100,1);

b_errall = b_res(:,3) + berr_res(:,3);
b_errall(:,2) =  b_res(:,3) - berr_res(:,3);
a_errall = [a_res(:,3) a_res(:,3)];


for i= 1:numel(X)
    yall = a_errall .* X(i).^(-b_errall +1);
    Y_errtop(i) = max(yall(:));
    Y_errbottom(i) = min(yall(:));
end
X = X';

% Load Josy Data
rockfalld = dlmread('data/rockfalls_Josy.dat');
VJ = sortrows(rockfalld);
timeJ = 1.5;
    
neventJ = numel(VJ);
numbersJ = [neventJ:-1:1]'/timeJ;

%load('results/Josy_fits.mat')
%Take ml values as reported in paper:
aJ_best = 47.9/timeJ; bJ_best = 1.67;
aJ_err = 1.3/timeJ; bJ_err = 0.07;
aJ_min = aJ_best+aJ_err; bJ_min = bJ_best-bJ_err;
aJ_max = aJ_best-aJ_err; bJ_max = bJ_best+bJ_err;

XJ = logspace( -2, 4, 100);
YJ_best = aJ_best * XJ.^(-bJ_best+1);
YJ_max = aJ_max * XJ.^(-bJ_max+1);
YJ_min = aJ_min * XJ.^(-bJ_min+1);

%% Check PDF notation

a_pdf = a_errall(:,1) ./ (b_res(:,3)-1);
vmin = 0.35; vmax = 130;
check_volume = @(a_pdf, b_pdf) a_pdf./(2 - b_pdf) .* (vmax.^(2-b_pdf) - vmin.^(2-b_pdf));
check_volume2 = @(a_pdf, b_pdf,vmin,vmax) a_pdf./(2 - b_pdf) .* (vmax.^(2-b_pdf) - vmin.^(2-b_pdf));
test_vol = @(data, vmin, vmax) sum(data(data>vmin & data < vmax));

checkvol = check_volume(a_pdf, b_res(:,3));

%% Plot Summary figure

figure('Name', 'Rockfall Events and Fit',...
    'Units', 'centimeters',...
    'Position', [1,1,20,14])



hold on
%plot this results with extended error range
[xpatch, ypatch] = create_patch(XJ, YJ_max, YJ_min);
patch('XData', xpatch, 'YData', ypatch,...
    'FaceColor',lightgreen,'FaceAlpha', 0.5,...
    'EdgeColor',darkgreen,'LineWidth',1.5,'LineStyle','-');
plot(XJ,YJ_best, 'LineWidth',2.5,'LineStyle','-','Color',darkgreen)

% [xpatch, ypatch] = create_patch(X, Y_top, Y_bottom);
% patch('XData', xpatch, 'YData', ypatch,...
%     'FaceColor',lightblue,'FaceAlpha', 0.3,...
%     'EdgeColor',darkblue,'LineWidth',2.0,'LineStyle','--');


[xpatch, ypatch] = create_patch(X, Y_errtop, Y_errbottom);
patch('XData', xpatch, 'YData', ypatch,...
    'FaceColor',lightblue,'FaceAlpha', 0.5,...
    'EdgeColor',darkblue,'LineWidth',2.5,'LineStyle','-');




%plot Events
mask_josy = VJ>=0.35;
loglog(VJ(mask_josy), numbersJ(mask_josy), 'd', 'MarkerSize', 6,...
    'Color', 0.5*[1,1,1], 'MarkerFaceColor',0.5*[1,1,1])
hold on
loglog(VJ(~mask_josy), numbersJ(~mask_josy), 'x', 'MarkerSize', 6,...
    'Color', 0.5*[1,1,1], 'MarkerFaceColor',0.5*[1,1,1])

mask_ro = (V >= V_ro1);
loglog(V(mask_ro), numbers(mask_ro), 'ko', 'MarkerSize', 6,...
    'MarkerFaceColor','k')
hold on
loglog(V(~mask_ro), numbers(~mask_ro), 'kx', 'MarkerSize', 6)


axis([0.05,1000,0.1,100])
box on

axmain = gca;
set (axmain, 'XScale', 'log')
set(axmain,'XTick', [0.01 0.1 1 10 100 1000])
set(axmain,'XTickLabel',{'0.01'; '0.1'; '1'; '10'; '100'; '1000'})
set (axmain, 'YScale', 'log')
set(axmain,'YTick',[0.1 1 10 100 ])
set(axmain,'YTickLabel',{'0.1'; '1'; '10';'100'})
axmain.LineWidth = 2.0;
axmain.FontSize = 16;
xlabel('rockfall volume [m^{3}]','FontWeight','bold','FontSize',18);
ylabel('cumulative number of events per year','FontWeight','bold','FontSize',18);

%legend
legh = legend({'Fit range by Strunden et al., 2015', ...
        'Fit range in this study', ...
        '1.5 year dataset', 'below rollover',...
        '5.2 year datset', 'below rollover'},'Location','southwest');
set(legh,'FontSize',15)

text(800, 55,sprintf('this study: N(V) = [%.1f - %.1f] V^{ -[%.2f - %.2f]}',...
    min(a_errall(:)), max(a_errall(:)),min(b_errall(:))-1,max(b_errall(:))-1),...
    'Color',darkblue,'FontSize',14,...
    'VerticalAlignment','top','HorizontalAlignment','right')
text(800, 85,sprintf('Strunden et. al: N(V) = (%.1f +- %.1f) V^{ -(%.2f +- %.2f)}',...
    aJ_best, aJ_err, bJ_best-1, bJ_err),...
    'Color',darkgreen,'FontSize',14,...
    'VerticalAlignment','top','HorizontalAlignment','right')

%print('Figures/Fit_Comparison', '-dpdf', '-r300')


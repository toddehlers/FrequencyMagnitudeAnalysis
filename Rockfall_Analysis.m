%% Rockfall Data Analysis Update April 2018

%add matlab routines provided by A. Clauset and the data director
addpath('routines')

rockfalldata = dlmread('data/rockfall_data_all.txt');
%rockfalldata = dlmread('data/rockfalls_Josy.dat');
rockfalldata = sortrows(rockfalldata);
time = 5.2;

% 3 Methods: NL-Fit the CCDF values both with and without error weights,
% and once following the Methods  following Clauset et al, SIAM, 2009

V = rockfalldata(:,1);
dV = rockfalldata(:,2);
%Pad rockfall errors to 5% if error is 0
dV(dV<1e-4) = V(dV<1e-4) * 0.05;


nevent = numel(V);
numbers = [nevent:-1:1]'/time;

%% NL-Fit, find rollover both ways


% find rollover, need to have at least 10 and not more than 40% of entire 
% dataset of rockfalls below rollover
min_event_ro = 0;
max_event_ro = round(nevent*0.85);

%start known from previous fits
st = [ 20; 0.8 ];
st_ro = [200, 0.15];

gof_arr = zeros(max_event_ro,2);
gof_arr_vals = {'adjrsquare', 'BIC'};


for i=1+min_event_ro:max_event_ro
    
    %goodness of fit by nl fit
    ft = fittype( 'power1' );
    opts_NL = fitoptions( 'Method', 'NonlinearLeastSquares', 'Display', 'notify', ...
            'StartPoint', st, 'Exclude', V<V(i));

    [fitres_bare, gof_bare] = fit(V, numbers, ft, opts_NL );
    
    %likelihood ratio after Clauset et al, 
    [~, ~, L, KS_gof] = plfit(V,'range',1.5:0.002:2.0,'xmin',V(i)); 
    
    
    gof_arr(i,1) = gof_bare.adjrsquare;
    gof_arr(i,2) = KS_gof;
    gof_arr(i,3) = L - 0.5*V(i)*(nevent-i);
end







%% Position of maxima in R2

%second peak in R2
[r2max2, i_r2max2] = max(gof_arr(:,1));
%first peak in R2 is lower
mask = V(1:max_event_ro) < 0.7;
[r2max1, i_r2max1] = max(gof_arr(mask,1));
%minimum between peaks, index needs to be adjusted because of slicing
[r2min1,i_r2min1] = min(gof_arr(i_r2max1:i_r2max2,1)); i_r2min1 = i_r2min1 + i_r2max1;

V_ro1 = V(i_r2max1);
V_ro2 = V(i_r2max2);

%% manual override for V_ro1 because of wide peak plateau
V_ro1 = 0.29; 
[~, i_ro1] = min(abs(V-V_ro1));

%% select Rollover after Clauset et al, 2009, using Bayesian Information Criterion (BIC)

%calculate vmin and b from data
[b, V_roKS, ~] = plfit(V,'range',1.5:0.05:2.0);
[~, i_roKS] = min(abs(V-V_roKS));

%% Plot adjusted R^2-value of fit as fuction of removing first events

figure('Name', 'Finding Rollover',...
    'Units', 'centimeters',...
    'Position', [1,1,18,14])

yyaxis left
plot(V(1:max_event_ro),gof_arr(:,1),'LineWidth',2.5)
line([V_ro1 V_ro1], [0.8,1.1],'LineStyle','--','Color',0.3*[1 1 1],'LineWidth',2.0)
text(V_ro1+0.05, 0.953, 'peak in R^{2} value', 'Rotation', 90, ...
    'FontSize', 15, 'Color', 0.3*[1 1 1])


ylabel('adjusted R^2 of nonlinear fit')
xlabel( 'Rollover Volume [m^3]')

yyaxis right
plot(V(1:max_event_ro),gof_arr(:,2),'LineWidth',2.5)
line([V_ro2 V_ro2], [0,1],'LineStyle','--','Color',0.3*[1 1 1],'LineWidth',2.0)
text(V_ro2-0.05, 0.055, 'lowest KS value', 'Rotation', 90, ...
    'FontSize', 15, 'Color', 0.3*[1 1 1])
ylabel('KS statistic')

ax=gca;

ax.YAxis(1).Limits = [0.9,1.01];
ax.YAxis(2).Limits =[0.025,0.15];
ax.XLim = [0,1.5];
set(ax,'XTick', 0:0.25:1.5)

ax.LineWidth = 2.0;
ax.FontSize = 16;

%print('Figures/Rollover_Limits', '-dpdf', '-r300')

%% create Range of possible rollovers, log spaced

nstep_ro = 25;
ro_range = exp(linspace(log(min(V_ro1,V_ro2)),log(max(V_ro1,V_ro2)),nstep_ro));

v_eval = 1000; %volume of event for return time calculation
t_eval = 100; %time interval for which the maximum recurring volume is calculate (e.g. 100-year-event)
vmin_er = V_ro1;  %lower volume limit for erosion integration
vmax_er = 1e4; %upper volume limit for erosion integration

%% Nonlinear Fit


%inital values for NL fit, known from previous fits
st = [ 20; 0.6 ];

a_res = zeros(nstep_ro,3);
aerr_res = zeros(nstep_ro,3);
b_res = zeros(nstep_ro,3);
berr_res = zeros(nstep_ro,3);
tret_res = zeros(nstep_ro,3);
erosion_res = zeros(nstep_ro,3);

%%
for i=1:nstep_ro
    % Volume and position of rollover
    V_ro = ro_range(i);
    [~, i_vro] = min(abs(V-V_ro));
    v_next = V(i_vro);
    n_vro = numbers(i_vro);

    %Linear fit of Log data
    [fit_log, gof_log] = fit(log(V),log(numbers),'poly1','Exclude',V<V_ro);

    %save in artificial fit object for evaluation of tret and erosion
    fit_linlog = cfit(fittype('power1'), exp(fit_log.p2), fit_log.p1);
    % caluctate t_return and average erosion
    [tret_log, ~, erosion_log] = eval_fit(fit_linlog, v_eval, t_eval, vmin_er, vmax_er);

    %save results in global array
    fitconf = confint(fit_log);
    a_res(i,1) = exp(fit_log.p2);
    aerr_res(i,1) = 0.5*(abs(diff(fitconf(:,2))));
    b_res(i,1) = 1-fit_log.p1;
    berr_res(i,1) = 0.5*(abs(diff(fitconf(:,1))));
    erosion_res(i,1) = erosion_log;
    tret_res(i,1) = tret_log;

    % The power law fit on the original data
    opts_NL = fitoptions( 'Method', 'NonlinearLeastSquares', 'Display', 'notify', ...
            'StartPoint', st, 'Exclude', V<V_ro); 
    fit_NL = fit(V, numbers, 'power1', opts_NL);

    %caluctate t_return and average erosion
    [tret_NL, ~, erosion_NL] = eval_fit(fit_NL, v_eval, t_eval, vmin_er, vmax_er);

    %save results in global array
    fitconf = confint(fit_NL);
    a_res(i,2) = fit_NL.a;
    aerr_res(i,2) = 0.5*(abs(diff(fitconf(:,1)))); %fitconf(:,1) here, different model!
    b_res(i,2) = 1-fit_NL.b;
    berr_res(i,2) = 0.5*(abs(diff(fitconf(:,2)))); %fitconf(:,2) here, different model!
    erosion_res(i,2) = erosion_NL;
    tret_res(i,2) = tret_NL;

    % Calculate fit parameters following Clauset et al, 2009

    %calculate vmin and b from data
    [b_ml, ~, ~] = plfit(V,'range',1.5:0.005:2.0,'xmin',V_ro);
    %determine b error
    [berr_ml, ~, ~] = plvar(V,'xmin',V_ro,'silent');
    %caclulate p value and goodness of fit
    [p_ml,gof_ml]=plpva(V, V_ro, 'silent');


    % save in cfit object
    a_ml = n_vro * V_ro^(b_ml-1);
    a_next = n_vro * v_next^(b_ml-1);

    fit_ml = cfit(fittype('power1'), a_ml, -b_ml+1);
    [tret_ml, ~, erosion_ml] = eval_fit(fit_ml,  v_eval, t_eval, vmin_er, vmax_er);

    %save results in global array
    a_res(i,3) = a_ml;
    aerr_res(i,3) = abs(a_ml-a_next);
    b_res(i,3) = b_ml;
    berr_res(i,3) = berr_ml;
    erosion_res(i,3) = erosion_ml;
    tret_res(i,3) = tret_ml;
end
%% adjust aerr

aerr_mean = mean(aerr_res./a_res);
for i=1:nstep_ro
    aerr_res(i,:) = a_res(i,:) .* aerr_mean;
end     




%% Save results

save('results/5years_all', 'V_ro1', 'V_ro2', 'V_roKS', 'i_r2min1', 'gof_arr', ...
     'a_res', 'aerr_res', 'b_res', 'berr_res', 'erosion_res', 'tret_res',   ...
     'X', 'Y_top', 'Y_bottom')

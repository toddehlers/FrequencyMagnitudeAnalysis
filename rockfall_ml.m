function [ fitparams, total_vol, time_return, v_ro ] = rockfall_ml ( rockfalls, time, v_ro)
%Calculate Maximum likelihood fit to rockfall inventory
% 
%This function calculates power law exponent b to the cumulative numbers of 
%rockfalls (N = a*V^-b), using the Maximum Likelihood method. The parameter
%a is estimated from linear regression of logarithmically transformed data.
%
%Input parameters:
% - rockfalls:  a (n x 1) table with rockfall volumes or a (n x 2) table of
%     rockfall volumes and their error (errors will be inored in
%     calculation)
% - time: the time in years over which the inventory was calculated
% - v_ro: the rollover volume in the rockfall inventory. If v_ro<0, the
%     rollover volume will be calculated by the routine fine_rollover using
%     standard parameters
% 
%Output:
% - fitparams is a 1x5 array with fit parameters:
%    (a, a_error, b, b_error, R^2-value of fit)
%    Keep in mind that the errors in a correspond to the 95%-confidence
%    interval of the fit and not the uncertainty in the data inventory.
% - total_vol (vmin, vmax) is a function that calculates the total eroded
%    volume caused by rockfalls of sizes ranging from vmin to vmax.
% - time_return (vmin) is a function that calculates the average 
%    recurrence time of any one event larger than vmax.
% - v_ro is the rollover volume used in calculations.

% Created by Matthias Schmiddunser-Nettesheim, April 2015
% eMail: matthias.schmiddunser@uni-tuebingen.de or todd.ehlers@uni-tuebingen.de

%-------------------------------------------------------------------------%
%% groom rockfall data for use
%-------------------------------------------------------------------------%

%check size of rockfall inventory
dim_rockf = size(rockfalls);
has_verr = min(dim_rockf);

%put rockfalls in row-form
if (dim_rockf(1) < dim_rockf(2)) 
    rockfalls = rockfalls';
end

%ignore errors if input has them
if (has_verr == 2)
    % Rockfalls have error associated with them
    rockfalls = rockfalls(:,1);
end

%if no rollover volume is given, calculate
if (v_ro < 0)
    v_ro = find_rollover(rockfalls,10,0.5);
end
 
%-------------------------------------------------------------------------%
%% Prepare data for fitting
%-------------------------------------------------------------------------%
 
%Create ccdf table
ccdf = transpose(numel(rockfalls):-1:1);
%normalize for time
ccdf = ccdf/time;

% select only values above rollover
fitmask = (rockfalls>=v_ro);
rockf = rockfalls(fitmask);
 
%-------------------------------------------------------------------------%
%% Apply Maximum Likelyhood method
%-------------------------------------------------------------------------%

% b is calculated by the formula given in Aki, 1965
b = 1/(log(10)*(mean(log10(rockf)) - log10(v_ro)));
sig = b/sqrt(numel(rockf));

% do fit to the CCDF with restricted b to get a
fitopt = fitoptions('poly1');
fitopt.Lower = [-b, -100];
fitopt.Upper = [-b, 100];
fitx = log(rockf);
fity = log(ccdf(fitmask));

[linfit, lingof] = fit(fitx,fity,'poly1',fitopt);
fitconf = confint(linfit);

%best result parameters for form N = a * vol ^ -b
a = exp(linfit.p2);

%Create Output array
fitparams = [ a, 0.5*abs(exp(fitconf(2,2))-exp(fitconf(1,2))), ...
               b, sig, lingof.rsquare ];
        
%-------------------------------------------------------------------------%
%% Defining evaluation functions
%-------------------------------------------------------------------------%

%return time is the inverse of number of events per year for a given volume
%range
time_return = @ (volmin)  volmin^(b)/a;

%Total eroded volume is integral over probability density times volume
total_vol = @ (volmin, volmax) a*b/(1-b) * (volmax^(1-b) - volmin^(1-b));

end

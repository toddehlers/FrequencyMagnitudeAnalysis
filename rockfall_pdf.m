function [ fitparams, total_vol, time_return, v_ro, binout ] = rockfall_pdf ( rockfalls, time, v_ro, nbins )
%Calculate PDF fit to rockfall inventory
%
%This function calculates the frequency density (dN/dV) of rockfalls, 
%creates a histogram of the logarithmically transformed frequency density 
%and fits a power law (dN/dV = a*V^-b) to the histogram.
%
%Input parameters:
% - rockfalls:  a (n x 1) table with rockfall volumes or a (n x 2) table of
%     rockfall volumes and their error
% - time: the time in years over which the inventory was calculated
% - v_ro: the rollover volume in the rockfall inventory. If v_ro<0, the
%     rollover volume will be calculated by the routine fine_rollover using
%     standard parameters
% - nbins: number of bins in which the rockfalls between v_ro and v_max
%     will be sorted. Output has additional bins for events under v_ro
% 
%Output:
% - fitparams is a 1x5 array with fit parameters:
%    (a, a_error, b, b_error, R^2-value of fit)
%    Keep in mind that the errors correspond to the 95%-confidence interval
%    of the fit and not the uncertainty in the data inventory.
% - total_vol (vmin, vmax) is a function that calculates the total eroded
%    volume caused by rockfalls of sizes ranging from vmin to vmax
% - time_returs (vmin, vmax) is a function that calculates the average 
%    recurrence time of any one event ranging from vmin to vmax
% v_ro is the rollover volume used in calculations
% binout is the histogram used for fitting:
%    -1st column: mean volume of bin
%    -2nd column: event frequency density per year
%    -3rd column: bin was used(1) for fitting or discarded (0)

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

%Sort rockfall inventory and determine maximum volume for binning
if (has_verr == 1)
    % Rockfalls have no error associated with them
    rockfalls = sort(rockfalls);
    vmax = rockfalls(end)*1.1;
    vmin = rockfalls(1)*0.9;
    
    %if no rollover volume is given, calculate
    if (v_ro < 0)
        v_ro = find_rollover(rockfalls,10,0.5);
    end
    
    %create sub-inventory of events less and larger than v_ro
    rockfex = rockfalls(rockfalls < v_ro);
    rockf = rockfalls(rockfalls >= v_ro);

elseif (has_verr == 2)
    % Rockfalls have error associated with them
    rockfalls = sortrows(rockfalls,1);
    vmax = rockfalls(end,1)+2*rockfalls(end,2); %2 sigma distance
    vmin = rockfalls(1,1)-2*rockfalls(1,2);
    
    %if no rollover volume is given, calculate
    if (v_ro < 0)
        v_ro = find_rollover(rockfalls(:,1),10,0.5);
    end
end

%-------------------------------------------------------------------------%
%% Binninng the data
%-------------------------------------------------------------------------%

%determine bin width (in log space)
wbin = (log(vmax)-log(v_ro))/nbins;
bins = zeros(nbins,1);

%bin boundaries in log space
binedge = log(v_ro)+ (0:nbins)'*wbin;
%bin position is mean in log space
binpos = exp(0.5*(binedge(2:end)+binedge(1:end-1)));

if (has_verr == 1)
    % Count events
    binid = ceil( (log(rockf)-log(v_ro)) / wbin );
    for i=1:nbins
        bins(i) = sum(binid==i);
    end
elseif (has_verr == 2)
     % Events are normal distribution functions rather than single points
     % bin count is determined by integrating over all events from lower
     % to upper bin boundary
     for i=1:nbins
         lbin = exp(binedge(i)); ubin = exp(binedge(i+1));
         for k=1:size(rockfalls,1);
             bins(i) = bins(i) + ...
                 normcdf(ubin,rockfalls(k,1),rockfalls(k,2)) - ...
                 normcdf(lbin,rockfalls(k,1),rockfalls(k,2));
             
         end
     end         
end

%normalize bins for time in years
bins = bins/time;

%normalize bins by bin width
PDF = bins./diff(exp(binedge));

%-------------------------------------------------------------------------%
%% Fitting the PDF
%-------------------------------------------------------------------------%


% Linear fit in Log-Space
% Needs to exclude empty bins
fitx = log(binpos(PDF>0));
fity = log(PDF(PDF>0));
[linfit, lingof] = fit(fitx,fity,'poly1');
fitconf = confint(linfit);

%best result parameters for form dN_R = a * vol ^ -b
a = exp(linfit.p2);
b = -linfit.p1;

%Create Output array
fitparams = [ a, 0.5*abs(exp(fitconf(2,2))-exp(fitconf(1,2))), ...
               b, 0.5*abs(fitconf(1,1)-fitconf(2,1)), lingof.rsquare ];
           
%-------------------------------------------------------------------------%
%% Defining evaluation functions
%-------------------------------------------------------------------------%
% See Barlow et al., 2012 for derivation of thess formulas or do it by hand
% on paper

%return time is the inverse of number of events per year for a given volume
%range
time_return = @ (volmin, volmax) 1/(a/(1-b) * (volmax^(1-b)-volmin^(1-b)));

%Total eroded volume is integral over probability density times volume
total_vol = @ (volmin, volmax) a/(2-b) * (volmax^(2-b) - volmin^(2-b));

%-------------------------------------------------------------------------%
%% Complete the binning for events below v_ro
%-------------------------------------------------------------------------%

%determine number of bins (in log space)
nbinex = ceil((log(v_ro)-log(vmin))/wbin);
vbinmin = v_ro*exp(-nbinex*wbin);
binex = zeros(nbinex,1);


%bin boundaries in log space
binedgeex = log(vbinmin)+ (0:nbinex)'*wbin;
%bin position is mean in log space
binposex = exp(0.5*(binedgeex(2:end)+binedgeex(1:end-1)));

if (has_verr == 1)
    % Count events
    binidex = ceil( (log(rockfex)-log(vbinmin)) / wbin );
    for i=1:nbinex
        binex(i) = sum(binidex==i);
    end
elseif (has_verr == 2)
     % Events are normal distribution functions rather than single points
     % bin count is determined by integrating over all events from lower
     % to upper bin boundary
     for i=1:nbinex
         lbin = exp(binedgeex(i)); ubin = exp(binedgeex(i+1));
         for k=1:size(rockfalls,1);
             binex(i) = binex(i) + ...
                 normcdf(ubin,rockfalls(k,1),rockfalls(k,2)) - ...
                 normcdf(lbin,rockfalls(k,1),rockfalls(k,2));
             
         end
     end         
end

%normalize bins for time in years
binex = binex/time;

%normalize bins by bin width
PDFex = binex./diff(exp(binedgeex));

fitusage = horzcat(zeros(1,nbinex), ones(1,nbins));

binout = transpose( [ binposex', binpos'; PDFex', PDF'; fitusage ] );

end

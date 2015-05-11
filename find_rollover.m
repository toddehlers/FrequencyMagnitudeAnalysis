function [ rollover ] = find_rollover ( rockfalls, skip, maxpos )
%Find the rollover in a rockfall inventory
%
% This function finds the rollover in a rockfall inventory by splitting it
% in two parts and fitting each subset. The split with the highest
% combined R^2-value to both fits is seen as best split.
%
% Input:
% - rockfalls:  a (n x 1) table with rockfall volumes or a (n x 2) table of
%     rockfall volumes and their error (errors will be inored in
%     calculation)
% - skip: the minimum number of events below the rollover
% - maxpos: the maximum fracion of events in the data repository below the 
%     the rollover volume. E.g., a value of 0.4 means that the rollover
%     must occur within the first 40% of recorded events.
% skip and maxpos depend on the data repository and must be chosen
% accordingly. Starting values could be skip=10 and maxpos=0.5.
%
% Created by Matthias Schmiddunser-Nettesheim, April 2015
% eMail: matthias.schmiddunser@uni-tuebingen.de or todd.ehlers@uni-tuebingen.de

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

% sort rock fall data
rockfalls = sort(rockfalls);

%preallocate
nevent = numel(rockfalls);
rockfalls = reshape(rockfalls, [nevent, 1]);
valuetable = zeros(round(nevent*maxpos)-skip, 5);

%crate ccdf
ccdf = transpose(nevent:-1:1);

%transform data into logspace
fitx = log(rockfalls);
fity = log(ccdf);

for i=1+skip:round(nevent*maxpos)
    %do fits on subsets
    [linfit1, lingof1] = fit(fitx(1:i),fity(1:i),'poly1'); 
    [linfit2, lingof2] = fit(fitx(i+1:end),fity(i+1:end),'poly1');
    
    %write result in valuetable
    valuetable(i-skip,1) = linfit1.p2; %a1
    valuetable(i-skip,2) = linfit1.p1; %b1
    valuetable(i-skip,3) = linfit2.p2; %a2
    valuetable(i-skip,4) = linfit2.p1; %b2
    %product of rsquared values to evaluate fit quality
    valuetable(i-skip,5) = lingof1.rsquare * lingof2.rsquare;
end

%get position of highest product of r2
[~,ropos] = max(valuetable(:,5));

%extract fit parameters for this split
a1 = valuetable(ropos,1);
b1 = valuetable(ropos,2);
a2 = valuetable(ropos,3);
b2 = valuetable(ropos,4);

%find rollover as intercept of the two fit lines
rollover = exp((a1-a2)/(b2-b1));

end
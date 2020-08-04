function [xpatch, ypatch]  = create_patch(x, y1, y2)
% creates a two concatenated arrays of min and max to fill a patch in a
% plot.

% if (numel(x1) /= numel(x2))
%     exit
% end

shape = size(x); 
if (shape(1) > shape(2))
    x = x';
end

shape = size(y1); 
if (shape(1) > shape(2))
    y1 = y1';
end

shape = size(y2); 
if (shape(1) > shape(2))
    y2 = y2';
end

maxval = max([y1;y2]);
minval = fliplr(min([y1;y2]));

xpatch = [x, fliplr(x)];
ypatch = [maxval, minval];


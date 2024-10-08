function [alpha] = kripAlpha(dat, scale, makeplot)
% Calculate Krippendorff's alpha using original approach.
% Use as:
%   alpha = kripAlpha(X, method, makeplot)
% Where 
%   dat:      N observers x M observations. For time series M = t.
%   scale:    The method for calculating the error (i.e. delta) for
%             Krippendorpf's Alpha. Can be NOMINAL, ORDINAL, INTERVAL,
%             ANGLE or RATIO.
%   makeplot: Plot the data, 1 or 0 (default = 0).

% Calculate alpha with hist3 approach (absolute values)

% Check inputs
if nargin < 2
    error('Error: must have argument SCALE')
elseif nargin < 3
    makeplot = 0;
end

fprintf('This dataset has %i observers and %i data points.\n',size(dat, 1) , size(dat, 2) )
scale = lower(scale);

if ~any(contains({'nominal','ordinal','interval','angle','ratio'}, scale))
    error('Unknown scale of measurement');
end

% Get variables
allvals = unique(dat(~isnan(dat)));
if isa(allvals, 'logical'); allvals = int8(allvals); end
tdim = size(dat, 2);

Y = repmat(1:size(dat,2), size(dat,1), 1);
XY = [dat(:), Y(:)];
dO = hist3(XY, {allvals, 1:tdim});
dE = sum(dO,2);
n__ = sum(dO(:));              %length(dat(:));
nu_ = sum(dO,1);                 %size(dat,1);

% Optional plot
if makeplot
    figure;
    subplot(1,20,1:3); plot(dE, allvals); axis tight; hold on
    subplot(1,20,5:20); imagesc(1:size(dat,2),allvals,dO); colorbar
    set(gca,'YDir','normal')
end

% nominator
Zu = zeros(1, tdim);
for tt = 1:tdim
    if ~(nu_(tt) > 1)
        continue
    end
    % Find values
    Znucnuk = zeros(1, length(allvals)-1);  % Initiate with zeroes.
    vidx = find(dO(:,tt));                  % Find only index with observed values, all else will be zero anyway
    
    for ii = 1:length(vidx)-1
        idx = vidx(ii);
        c = allvals(idx);
        kvals = allvals(vidx(ii+1:end));

        switch scale
            case 'nominal'
                deltas = delta_nominal(c, kvals);
            case 'ordinal'
                deltas = delta_ordinal(dE, vidx, ii, kvals);               
            case 'interval'
                deltas = delta_interval(c, kvals);
            case 'angle'
                deltas = delta_angle(c, kvals);
            case 'ratio'
                deltas = delta_ratio(c, kvals);
        end
        nuc = dO(idx,tt);
        nuk = dO(vidx(ii+1:end), tt);
        Znucnuk(idx) = sum(nuc*nuk.*deltas);    
    end   
    Zu(tt) = sum(Znucnuk./(nu_(tt)-1));
end

% denominator
Zncnk = zeros(1, length(allvals));
for ii = 1:length(allvals)-1
    c = allvals(ii);
    kvals = allvals(ii+1:end);
    vidx = 1:length(allvals);

    switch scale
        case 'nominal'
            deltas = delta_nominal(c, kvals);
        case 'ordinal'
            deltas = delta_ordinal(dE, vidx, ii, kvals);               
        case 'interval'
            deltas = delta_interval(c, kvals);
        case 'angle'
            deltas = delta_angle(c, kvals);
        case 'ratio'
            deltas = delta_ratio(c, kvals);
    end        
    n_c = dE(ii);
    n_k = dE(ii+1:end);
    Zncnk(ii) = sum(n_c*n_k.*deltas);
end

Do = sum(Zu);
De = sum(Zncnk) / (n__-1);

alpha = 1 - (Do/De);

% Error functions
function deltas = delta_nominal(c, kvals)
    deltas = ~(kvals==c);
end
function deltas = delta_ordinal(dE, vidx, ii, kvals)
    deltas = nan(length(kvals),1);
    kidx = vidx(ii+1:end);
    for gg = 1:length(kidx)
        deltas(gg) = (sum(dE(vidx(ii):kidx(gg))) - (dE(vidx(ii))+dE(kidx(gg)))/2)^2;
    end  
end
function deltas = delta_interval(c, kvals)
    deltas = (kvals-c).^2;
end
    function deltas = delta_angle(c, kvals)
    deltas = sin((kvals-c)/2).^2;
end
function deltas = delta_ratio(c, kvals)
    deltas = ((c-kvals)./(c+kvals)).^2;
end

end
%END
function [alpha, cfg] = kripAlpha(dat, scale)
% Calculate Krippendorff's alpha using original approach.
% Use as:
%   [alpha, cfg] = kripAlpha(X, method)
% Where 
%   dat:      N observers x M observations. For time series M = t. Data
%             must be numeric.
%   scale:    The method for calculating the error (i.e. delta) for
%             Krippendorpf's Alpha. Can be NOMINAL, ORDINAL, INTERVAL,
%             ANGLE_DEG, ANLGE_RAD or RATIO.
%
% Output:
%   alpha:    The Alpha value
%   cfg:      Setting for bootstrap procedure.

% Calculate alpha with hist approach (absolute value)

% Check inputs
if nargin < 2
    error('Error: must have argument SCALE')
end

scale = lower(scale);

if ~any(contains({'nominal','ordinal','interval','angle','angle_rad','angle_deg','ratio'}, scale))
    error('Unknown scale of measurement');
end

if strcmp(scale, 'angle')
    scale = 'angle_deg';
end

% Get variables
fprintf('This dataset has %i observers and %i data points.\n',size(dat, 1) , size(dat, 2) )

allvals = unique(dat(~isnan(dat)));
if isa(allvals, 'logical'); allvals = int8(allvals); end

% Y = repmat(1:size(dat,2), size(dat,1), 1);
% XY = [dat(:), Y(:)];
% dO = hist3(XY, {allvals, 1:size(dat, 2)});
% dE = sum(dO(:,sum(dO)>1),2);

dE = hist(dat(:), allvals)';        % Expected count
% nu_ = sum(dO, 1);
nu_ = sum(~isnan(dat));
n__ = sum(nu_(nu_>1));
% clear XY Y

% f = waitbar(0,'Calculating...');

% nominator
Zu = 0;
for tt = 1:length(nu_)
    if ~(nu_(tt) > 1)
        continue
    end

    % Find values
    dOu = hist(dat(:,tt),allvals);
    vidx = find(dOu);

    Znucnuk = 0;            % Initiate as zero.
    deltas = zeros(length(vidx) -1, 1);   % Preallocate
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
            case 'angle_deg'
                deltas = delta_angle_deg(c, kvals);
            case 'angle_rad'
                deltas = delta_angle_rad(c, kvals);
            case 'ratio'
                deltas = delta_ratio(c, kvals);
        end
        nuc = dOu(idx); %sum(dat(:,tt) == c); % dO(idx,tt);
        nuk = dOu(vidx(ii+1:end))';
        Znucnuk = Znucnuk + sum(nuc*nuk.*deltas);
    end   

%     Zu(tt) = sum(Znucnuk./(nu_(tt)-1));
    Zu = Zu + sum(Znucnuk./(nu_(tt)-1));
%     dt = toc;
%     disp(dt)
end

% denominator
Zncnk = 0;
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
        case 'angle_deg'
            deltas = delta_angle_deg(c, kvals);
        case 'angle_rad'
            deltas = delta_angle_rad(c, kvals);
        case 'ratio'
            deltas = delta_ratio(c, kvals);
    end        
    n_c = dE(ii);
    n_k = dE(ii+1:end);
    Zncnk = Zncnk + sum(n_c*n_k.*deltas);
end

alpha = 1 - (Zu/Zncnk) * (n__-1);

% close(f)

% Variables for bootstrapping
cfg.n__     = n__;
cfg.mu      = nu_;
cfg.dE      = dE;
cfg.allvals = allvals;
cfg.scale   = scale;

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
function deltas = delta_angle_deg(c, kvals)
    deltas = sin(pi*(c-kvals)/360).^2;
end
function deltas = delta_angle_rad(c, kvals)
    deltas = sin(pi*(c-kvals)/(2*pi)).^2;
end
function deltas = delta_ratio(c, kvals)
    deltas = ((c-kvals)./(c+kvals)).^2;
end

end
%END
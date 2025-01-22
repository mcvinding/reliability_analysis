function [bootalphas] = bootstrap_alpha(dat, cfg, nboot)
% Calculate bootstrap distribution for Krippendorff's alpha.
% Use as:
%   bootalphas = bootstrap_alpha(dat, cfg, nboot)
% Where 
%   dat:      N observers x M observations. For time series M = t.
%   cfg:      Output from main functions.
%   nboot:    Number of bootstraps (default = 1000).

if nargin < 3
    nboot = 1000;
end

% Init values
n__    = cfg.n__;
mu     = cfg.mu;
vals   = cfg.allvals;
scale  = cfg.scale;
dE     = cfg.dE; 

No = sum(((mu-1).*mu)/2);

% Calculate traditional D_e
Zncnkd = 0;
for i = 1:length(dE)
    c = vals(i);
    nc = dE(i);

    switch scale
        case 'nominal'
            deltas = delta_nominal(c, vals);
        case 'ordinal'
            deltas = zeros(size(dE));
            for g = 1:length(dE)
                if g < i
                    deltas(g) = (sum(dE(g:i)) - (dE(i)+dE(g))/2)^2;
                else
                    deltas(g) = (sum(dE(i:g)) - (dE(i)+dE(g))/2)^2;
                end
            end
        case {'interval','n2fast','alphaprime','prime'}
            deltas = delta_interval(c, vals);
        case {'angle', 'angle_deg'}
            deltas = delta_angle_deg(c, vals);
        case 'angle_rad'
            deltas = delta_angle_rad(c, vals);
        case 'ratio'
            deltas = delta_ratio(c, vals);
    end
    Zncnkd = Zncnkd + sum(nc*dE.*deltas);
end

De = Zncnkd / (n__*(n__-1));

%% Get all pairs
% For N=2 this will be the same as actual data.
if strcmp(scale, 'n2fast')
    pairs = dat';
% elseif any(strcmp(scale, {'alphaprime','prime'}))
else
    pairs = [];
    for u = 1:length(mu)
        for i = 1:length(dat(:,u))-1
            for j = (i+1):length(dat(:,u))
                if ~(dat(i,u)==dat(j,u)) && ~(isnan(dat(i,u))) && ~(isnan(dat(j,u)))
                    pairs = [pairs; [dat(i,u), dat(j,u)]];
                end
            end
        end
    end
end
clear dat

%% Make a sample of deltas
delta2 = zeros(length(pairs),1);
for p = 1:length(pairs)
    switch scale
        case 'nominal'
            delta2(p) = delta_nominal(pairs(p,1), pairs(p,2));
        case 'ordinal'
            i = find(vals==pairs(p,1));
            g = find(vals==pairs(p,2));
            if g < i
                delta2(p) = (sum(dE(g:i)) - (dE(i)+dE(g))/2)^2;
            else
                delta2(p) = (sum(dE(i:g)) - (dE(i)+dE(g))/2)^2;
            end
        case {'interval','n2fast','alphaprime','prime'}
            delta2(p) = delta_interval(pairs(p,1), pairs(p,2));
        case {'angle', 'angle_deg'}
            delta2(p) = delta_angle_deg(pairs(p,1), pairs(p,2));
        case 'angle_rad'
            delta2(p) = delta_angle_rad(pairs(p,1), pairs(p,2));
        case 'ratio'
            delta2(p) = delta_ratio(pairs(p,1), pairs(p,2));
    end
end
clear pairs

%% Do bootstrapping
bootalphas = zeros(nboot,1);
for b = 1:nboot
    dsum = 0;
    for u = 1:length(mu)
        pairr = (mu(u)-1)*mu(u)/2;
        for ii = 1:pairr
            r = randi(No);                  % Random sample
            if r <= length(delta2)
                Er = 2 * delta2(r) / (n__ * De);
                d = Er/(mu(u)-1);
                dsum = dsum + d;
            end
        end
    end
    bootalphas(b) = 1-dsum;
end

% Error functions
function deltas = delta_nominal(c, kvals)
    deltas = ~(kvals==c);
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

function [alpha, boot] = reliability_analysis(X, method, bootstrap)
% Caluclate Krippendorf's Alpha.
% Use as:
%   [alpha, boot] = reliability_analysis(X, method, bootstrap)
% Where 
%   X:          N observers x M observations. For time series M = t.
%   method:     The method for calculating the error (i.e. delta) for
%               Krippendorpf's Alpha. Options can be:
%                 'interval': alpha for INTERVAL data.
%                 'ordinal': alpha for ORDINAL data.
%                 'nominal': alpha for NOMINAL data.
%                 'ratio' : alpha for RATIO data.
%                 'angle_deg' : alpha for PHASE data (in degrees).
%                 'angle_rad' : alpha for PHASE data (in radians).
%                 'n2fast_interval': faster computation of alpha for INTERVAL
%                   data with N=2 observers.
%                 'n2fast_nominal': faster computation of alpha for NOMINAL 
%                   data with N=2 observers.
%                 'alphaprime': approximation of the observation matrices 
%                   using calculation by binning. Suitable for large datasets 
%                   of INTERVAL data with arbitrary numerical precision.
%   bootstrap:  (int) Calculate bootstrapping distribution for calculating 
%               CI of alpha. Number of bootstrapping. 0 = do not calculate 
%               bootstrapping (default).

% This is a wrapper function, see KRIPALPHA, ALPHAPRIME, KRIPALPHAN2FAST,
%   and BOOTSTRAP_ALPHA for further documentation.

% Undocumented options:
%  method:      'prime': same as 'alphaprime'.
%               'n2fast': same as 'n2fast_interval'. Kept for backwards
%                 compatability.
%               'angle': same as 'angle_deg'.Kept for backwards compatability.

% Check input
if nargin < 2
    error('Error: must have argument SCALE')
elseif nargin < 3
    bootstrap = 0;
end

if bootstrap < 0
    error('Error: Number of bootstraps cannot be negative.')
end

% Check input data
if length(size(X)) > 2
    error('Error: input data should be a 2-dimensional NxM matrix. Input has %i dimensions.', length(size(X)))
end

% Check methos
method = lower(method);
if ~any(strcmp(method, {'alphaprime', 'prime', ...
        'interval', 'ordinal', 'nominal', 'angle', 'angle_deg', 'angle_rad', 'ratio', ...
        'n2fast', 'n2fast_interval', 'n2fast_nominal'}))
    error('Error: method \"%s\" is not supported', upper(method))
end

% Switch methods and run
switch(method)
    case {'alphaprime', 'prime'}
        disp('Calculating Alpha with denisty approximation...\n')
        [alpha, cfg] = alphaprime(X);
        
    case {'interval', 'ordinal', 'nominal', 'ratio', 'angle', 'angle_rad', 'angle_deg'}
        fprintf('Calculating Alpha for %s data with exact precision...\n', upper(method))
        [alpha, cfg] = kripAlpha(X, method);
        
    case {'n2fast', 'n2fast_interval', 'n2fast_nominal'}
        if strcmp(method, 'n2fast'); method = 'n2fast_interval'; end
        fprintf('Calculating Alpha with fast method for N=2\n')
        [alpha, cfg] = kripAlphaN2fast(X, method);
end

if bootstrap > 0
    fprintf('Running bootstrap procedure...\n')
    boot = bootstrap_alpha(X, cfg, bootstrap);

    ci = prctile(boot, [2.5, 97.5]); sig=0.8;
    fprintf('  Alpha = %.3f (95%%CI: %.3f-%.3f [%i bootstraps])\n', alpha, ci(1), ci(2), bootstrap)
    fprintf('  Probability of Alpha >= %.2f: %.3f (P-value = %.3f)\n', sig, mean(boot > 0.8), 1-mean(boot > 0.8));
end

%END
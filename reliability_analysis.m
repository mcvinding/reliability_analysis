function [alpha, boot] = reliability_analysis(X, method, bootstrap)
% Caluclate Krippendorf's Alpha.
% Use as:
%   [alpha, boot] = reliability_analysis(X, method, bootstrap)
% Where 
%   X:          N observers x M observations. For time series M = t.
%   method:     The method for calculating the error (i.e. delta) for
%               Krippendorpf's Alpha. Options can be:
%                 'interval': alpha for INTERVAL data using the exact method.
%                 'ordinal': alpha for ORDINAL data using the exact method.
%                 'nominal': alpha for NOMINAL data using the exact method.
%                 'angle' : alpha for PHASE data using the exact method.
%                 'ratio' : alpha for RATIO data using the exact method.
%                 'alphaprime': approximation of the observation matrices 
%                   using calculation by binning. Suitable for large datasets 
%                   of INTERVAL data with arbitrary numerical precision.
%                'n2fast': fast computation of alpha for INTERVAL data with
%                 N=2 observers and no missing cases.
%   bootstrap:  (int) Calculate bootstrapping distribution for calculating 
%               CI of alpha. Number of bootstrapping. 0 = do not calculate 
%               bootstrapping (default).

% This is a wrapper function, see KRIPALPHA, ALPHAPRIME, KRIPALPHAN2FAST
% for further documentation.

% Check input
if nargin < 2
    error('Error: must have argument SCALE')
elseif nargin < 3
    bootstrap = 0;
end

if bootstrap < 0
    error('Error: Number of bootstraps must be positive or zero for no bootstrapping.')
end

% Check input data
if length(size(X)) > 2
    error('Error: input data should be a 2-dimensional NxM matrix. Input has %i dimensions.', length(size(X)))
end

% Check methos
method = lower(method);
if ~any(strcmp(method, {'alphaprime', 'prime', 'interval', 'ordinal', 'nominal', 'n2fast', 'angle', 'ratio'}))
    error('Error: method \"%s\" is not supported', upper(method))
end

% Switch methods and run
switch(method)
    case {'alphaprime', 'prime'}
        disp('Calculating Alpha with denisty approximation...\n')
        [alpha, cfg] = alphaprime(X);
        
    case {'interval', 'ordinal', 'nominal', 'ratio', 'angle'}
        fprintf('Calculating Alpha for %s data with exact precision...\n', upper(method))
        [alpha, cfg] = kripAlpha(X, method);
        
    case 'n2fast'
        fprintf('Calculating Alpha with fast method for N=2 and interval data\n')
        [alpha, cfg] = kripAlphaN2fast(X);        
end

if bootstrap > 0
    fprintf('Running bootstrap procedure...\n')
    boot = bootstrap_alpha(X, cfg, bootstrap);
end

%END
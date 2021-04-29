function [alpha] = reliability_analysis(X, method)
% Caluclate Krippendorf's Alpha.
% Use as:
%   alpha = reliability_analysis(X, method)
% Where 
%   X:      N observers x M observations. For time series M = t.
%   method: The method for calculating the error (i.e. delta) for
%           Krippendorpf's Alpha. Options can be:
%            'interval': K's alpha for INTERVAL data using exact method.
%            'ordinal': K's alpha for ORDINAL data using exact method.
%            'nominal': K's alpha for NOMINAL data using exact method.
%            'alphaprime': approximation of coincidence matrices using
%               density calculation and interpolation. Only suitable for 
%               large datasets of INTERVAL data with arbitrary numerical 
%               precision.
%            'n2fast': fast computation of K's alpha for INTERVAL data with
%               N = 2 observers and no missing cases.
%

% This is a wrapper function, see KRIPALPHA, ALPHAPRIME, KRIPALPHAN2FAST
% for proper documentation.

% Check input data
if length(size(X)) > 2
    error('Error: input data should be a 2-dimensional NxM matrix. Input has %i dimensions.', length(size(X)))
end

% Check methos
method = lower(method);
if ~any(strcmp(method, {'alphaprime', 'interval', 'ordinal', 'nominal', 'n2fast'}))
    error('Error: method \"%s\" is not supported', upper(method))
end

% Switch methods and run
switch(method)
    case 'alphaprime'
        disp('Calculating Alpha with denisty approximation...')
        alpha = alphaPrime(X);
        
    case {'interval', 'ordinal', 'nominal'}
        fprintf('Calculating K´s Alpha for %s data with exact precision...\n', upper(method))
        alpha = kripAlpha(X, method);
        
    case 'n2fast'
        fprintf('Calculating K´s Alpha with fast method for N=2 and interval data')
        alpha = kripAlphaN2fast(X);
end
        
%END
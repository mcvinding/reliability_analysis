function [alpha] = reliability_analysis(X, method)
% Use as:
%   alpha = reliability_analysis(X, method)
% Where 
%   X:      N observers x M observations. For time series M = t
%   method: The method for calculating the error (i.e. delta) for
%           Krippendorpf's Alpha. Can be ... or 'AlphaPrime' for density
%           approximation (suitable for large datasets of interval data with
%           arbitrary numerical precision).


% Check input data
if length(size(X)) > 2
    error('Error: input data should be a 2-dimensional NxM matrix. Input has %i dimensions.', length(size(X)))
end

% Check methos
method = lower(method);
if ~any(strcmp(method, {'alphaprime', 'interval', 'ordinal', 'nominal'}))
    error('Error: method \"%s\" is not supported', upper(method))
end

% Switch methods and run
switch(method)
    case 'alphaprime'
        disp('Calculating Alpha with denisty approximation...')
        alpha = alphaPrime(X);
        
    case {'interval', 'ordinal', 'nominal'}
        fprintf('Calculating K´s Alpha for %s data with finite precision...\n', upper(method))
        alpha = kripAlpha(X, method);
end
        

end)
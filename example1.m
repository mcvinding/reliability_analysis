% Use example data frame from Krippendorff

A = [1 2 3 3 2 1 4 1 2 nan nan nan];
B = [1 2 3 3 2 2 4 1 2 5 nan 3];
C = [nan 3 3 3 2 3 4 2 2 5 1 nan];
D = [1 2 3 3 2 4 4 1 2 5 1 nan];

dat = [A; B; C; D];

%% Get Krippendorff's Alpha
alpha_int = kripAlpha(dat, 'interval'); disp('done')  % Assuming interval data
alpha_ord = kripAlpha(dat, 'ordinal'); disp('done')   % Assuming ordinal data

%% Approximation method (only interval and large datasets)
% Does not really make sense for this data.
alphap = alphaprime(dat); disp('done')                % Assuming interval data

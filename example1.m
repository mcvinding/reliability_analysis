% Use example dataframe from Krippendorff (ref)

A = [1 2 3 3 2 1 4 1 2 nan nan nan];
B = [1 2 3 3 2 2 4 1 2 5 nan 3];
C = [nan 3 3 3 2 3 4 2 2 5 1 nan];
D = [1 2 3 3 2 4 4 1 2 5 1 nan];

dat = [A; B; C; D];

%% Get Krippendorff's Alpha
alpha_int = kripAlpha(dat, 'interval');  % Assuming interval data
alpha_ord = kripAlpha(dat, 'ordinal');   % Assuming ordinal data

%% Approximation method
% Only interval and large datasets. Does not really make sense for this data.
alphap = alphaprime(dat);                % Assuming interval data

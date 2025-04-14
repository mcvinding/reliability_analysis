% Example dataframe from Krippendorff (2011)

A = [1 2 3 3 2 1 4 1 2 nan nan nan];
B = [1 2 3 3 2 2 4 1 2 5 nan 3];
C = [nan 3 3 3 2 3 4 2 2 5 1 nan];
D = [1 2 3 3 2 4 4 1 2 5 1 nan];

dat = [A; B; C; D];

%% Get Krippendorff's Alpha
alpha_int = reliability_analysis(dat, 'interval');      % Assuming interval data
alpha_ord = reliability_analysis(dat, 'ordinal');       % Assuming ordinal data
alpha_nom = reliability_analysis(dat, 'nominal');       % Assuming nominal data
alpha_rat = reliability_analysis(dat, 'ratio');         % Assuming ratio data
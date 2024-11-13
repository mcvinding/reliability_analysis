%% Example #2: Dataframe from Hayes & Krippendorff (2007)
% Get alpha assuming data is INTERVAL, ORDINAL, or NOMINAL. Use bootstrap 
% method to get 95% CI for alpha. Dataframe from Hayes & Krippendorff (2007)

A = [1 1 2 nan 0 0 1 1 2 2 nan 0 1 3 1 1 2 1 1 0 0 0 2 0 1 0 0 1 1 1 1 2 2 3 2 2 2 2 2 1];
B = [1 1 3 0 0 0 0 nan 2 1 1 0 2 3 1 1 1 2 1 0 0 0 3 0 2 1 0 2 1 1 1 1 2 2 2 2 2 2 2 1];
C = [2 0 3 0 0 0 2 2 2 1 0 0 2 2 1 1 2 3 0 0 1 nan 3 0 nan 1 0 1 2 2 0 2 nan 2 2 3 2 nan 2 1];
D = [nan 1 3 nan nan nan nan 0 nan 1 0 0 2 2 nan nan nan 3 1 nan 1 0 3 0 2 1 1 2 2 nan nan 1 2 2 nan nan nan 1 2 nan];
E = [2 nan nan 0 0 0 1 nan 2 nan nan nan nan 3 1 1 2 nan nan 0 nan 0 nan nan 2 nan 0 nan nan 2 0 nan 2 nan 2 2 2 2 nan 1];

dat = [A; B; C; D; E];

%% Get Krippendorff's Alpha

[alpha_int, boot_int] = reliability_analysis(dat, 'interval', 20000);  % Assuming interval data
[alpha_ord, boot_ord] = reliability_analysis(dat, 'ordinal', 20000);   % Assuming ordinal data
[alpha_nom, boot_nom] = reliability_analysis(dat, 'nominal', 20000);   % Assuming nominal data

%% Summary and plot
sig = 0.8;                          % Critical cutoff for "significance"

% Interval data
disp('INTERVAL:')
ci = prctile(boot_int, [2.5, 97.5]);
fprintf('Alpha = %.3f (CI: %.3f-%.3f)\n', alpha_int, ci(1), ci(2))
fprintf(['Probability of alpha being above threshold of %.2f:\n' ...
         '      P = %.3f\n'], sig, mean(boot_int < 0.8));
figure; histogram(boot_int, 30, 'Normalization', 'pdf');
xline(alpha_int, 'k', 'LineWidth',2);
xline(sig, 'r--', 'LineWidth',2);

% Ordinal data
disp('ORDINAL:')
ci = prctile(boot_ord, [2.5, 97.5]);
fprintf('Alpha = %.3f (CI: %.3f-%.3f)\n', alpha_ord, ci(1), ci(2))
fprintf(['Probability of alpha being above threshold of %.2f:\n' ...
         '      P = %.3f\n'], sig, mean(boot_ord < 0.8));
figure; histogram(boot_ord, 30, 'Normalization', 'pdf');
xline(alpha_ord, 'k', 'LineWidth',2);
xline(sig, 'r--', 'LineWidth',2);

% Nominal data
disp('NOMINAL:')
ci = prctile(boot_nom, [2.5, 97.5]);
fprintf('Alpha = %.3f (CI: %.3f-%.3f)\n', alpha_nom, ci(1), ci(2))
fprintf(['Probability of alpha being above threshold of %.2f:\n' ...
         '      P = %.3f\n'], sig, mean(boot_nom < 0.8));
figure; histogram(boot_nom, 30, 'Normalization', 'pdf');
xline(alpha_nom, 'k', 'LineWidth',2);
xline(sig, 'r--', 'LineWidth',2);

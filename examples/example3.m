%% Example #3: Test agreement between two time-series: Sine-waves with random noise. 
close all

%% Generate two time-series
StopTime    = 2;                 % seconds 
fs          = 1000;              % Sampling frequency (samples per second) 
err         = 0.35;               % Noise scaling
dt          = 1/fs;              % seconds per sample 
t           = (dt:dt:StopTime);  % seconds 
Freq        = 3;                 % Sine wave frequency (Hz) 

x = sin(2*pi*Freq*t)+randn(1,length(t))*err;
y = sin(2*pi*Freq*t)+randn(1,length(t))*err;
dat = [x; y];   % N observers vs M samples

plot(dat')

%% Get Alpha
[alpha_int, boot_int] = reliability_analysis(dat, 'interval', 20000);   % Exact method
[alpha_n2f, boot_n2f] = reliability_analysis(dat, 'n2fast', 20000);     % N=2 method (same a exact in this case)
[alpha_pri, boot_pri] = reliability_analysis(dat, 'alphaprime', 20000); % Alphaprime calculation with binning

%% Summary and plot
sig = 0.8;                          % Critical cutoff for "significance"

% Exact method
disp('EXACT METHOD:')
ci = prctile(boot_int, [2.5, 97.5]);
fprintf('Alpha = %.3f (CI: %.3f-%.3f)\n', alpha_int, ci(1), ci(2))
fprintf(['Probability of alpha being above threshold of %.2f:\n' ...
         '      P = %.3f\n'], sig, mean(boot_int < 0.8));
figure; histogram(boot_int, 30, 'Normalization', 'pdf');
xline(alpha_int, 'k', 'LineWidth',2);
xline(sig, 'r--', 'LineWidth',2);

% N=2 method
disp('N=2 METHOD:')
ci = prctile(boot_n2f, [2.5, 97.5]);
fprintf('Alpha = %.3f (CI: %.3f-%.3f)\n', alpha_n2f, ci(1), ci(2))
fprintf(['Probability of alpha being above threshold of %.2f:\n' ...
         '      P = %.3f\n'], sig, mean(boot_n2f < 0.8));
figure; histogram(boot_n2f, 30, 'Normalization', 'pdf');
xline(alpha_n2f, 'k', 'LineWidth',2);
xline(sig, 'r--', 'LineWidth',2);

% Alpha prime
disp('ALPHA PRIME METHOD:')
ci = prctile(boot_pri, [2.5, 97.5]);
fprintf('Alpha = %.3f (CI: %.3f-%.3f)\n', alpha_pri, ci(1), ci(2))
fprintf(['Probability of alpha being above threshold of %.2f:\n' ...
         '      P = %.3f\n'], sig, mean(boot_pri < 0.8));
figure; histogram(boot_pri, 30, 'Normalization', 'pdf');
xline(alpha_pri, 'k', 'LineWidth',2);
xline(sig, 'r--', 'LineWidth',2);

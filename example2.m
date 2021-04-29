% Test agreement  between two time-series: Sine-waves with random noise

%% Generate two time-series : settings
fs          = 200;               % Sampling frequency (samples per second) 
dt          = 1/fs;              % seconds per sample 
StopTime    = 1;                 % seconds 
t           = (dt:dt:StopTime);  % seconds 
Freq        = 5;                 % Sine wave frequency (hertz) 

%% Test
for ii = 1:10
    fprintf('Run %i of 10... ', ii) 
    for aa = 1:10
        % Generate two time-series
        x = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;
        y = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;
        dat = [x; y];   % N observers vs M samples

        % Plot data
        % figure(1);
        subplot(2,5,aa); 
        plot(t, x); hold on
        % subplot(2,5,aa); 
        plot(t, y);
        title(aa)

        % Test with kriAplha function
        tic
        alpha1(aa,ii) = kriAlpha(dat, 'interval'); disp('done')
        dT1(aa,ii) = toc;

        % New implementation
        tic
        alpha2(aa,ii) = kripAlpha(dat, 'interval'); disp('done')
        dT2(aa,ii) = toc;

        % Approximation method
        tic
        alphap(aa,ii) = alphaprime(dat); disp('done')
        dTp(aa,ii) = toc;
    end
    disp('done')
end
disp('DONE')
%%
[alpha1; alpha2; alphap]

figure;
plot(1:10, alpha1, 'o-'); hold on
plot(1:10, alpha2, 'o-'); hold on
plot(1:10, alphap, 'o-');



mean([dT1; dT2; dTp], 2)
std([dT1; dT2; dTp]')




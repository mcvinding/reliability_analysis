% Test agreement  between two time-series: Sine-waves with random noise

% Generate two time-series
fs      = 200; % Sampling frequency (samples per second) 
dt      = 1/fs; % seconds per sample 
StopTime = 1; % seconds 
t = (dt:dt:StopTime); % seconds 
Freq = 5; % Sine wave frequency (hertz) 


for aa = 0:10

x = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;
y = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;

figure(1);
% subplot(2,5,aa); 
plot(t, x); hold on
% subplot(2,5,aa); 
plot(t, y);
title(aa)

dat = [x; y];   % N observers vs M samples

%% Test with kriAplha function
% tic
% alpha1(aa) = kriAlpha(dat, 'interval');
% toc
%% New method
% TO DO:
% * Option to specify resolution
allvals = unique(dat(:));
N = 100;        % 100 = Default in linspace




end
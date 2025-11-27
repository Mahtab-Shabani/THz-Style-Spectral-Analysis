clear all; clc; close all; 

load data.mat;
% EMG1=x1(1001:6000);
EMG=x1(1001:6000);
y=EMG;
d=dalt;
t=thalt;
a=alpha;
b=belta;
fs=256;
nfft=1024;
noverlap = 750;

%% Time-frequency feature extraction

signals = {d, t, a, b, EMG};
names   = {'Delta','Theta','Alpha','Beta','EMG'};

window = 256;
noverlap = 200;
nfft = 512;

figure;

for i = 1:length(signals)
    subplot(5,2,2*i-1)
    spectrogram(signals{i}, window, noverlap, nfft, fs);
    title([names{i} ' - Spectrogram']); colorbar;

    % Compute average T-F energy
    [S,F,T] = spectrogram(signals{i}, window, noverlap, nfft, fs);
    TF_energy = mean(abs(S),2);

    subplot(5,2,2*i)
    plot(F, TF_energy);
    xlabel('Hz'); ylabel('Energy');
    title([names{i} ' - Frequency Energy']);
    grid on;
end

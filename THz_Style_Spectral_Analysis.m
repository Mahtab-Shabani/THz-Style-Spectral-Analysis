% THz_Style_Spectral_Analysis.m
% Time-Frequency and THz-style spectral biomarker extraction
% Compatible with MATLAB R2014a
%
% What it does (summary):
%  - Loads data.mat (expects variables: dalt, thalt, alpha, belta, x1, x2, x12120x2Eedf)
%  - For each signal: computes spectrogram (STFT), mean spectrum (via pwelch),
%    extracts spectral biomarkers (entropy, centroid, bandwidth, flatness, peak freq, powers),
%    smooths mean spectrum (for FDA-like handling), computes derivative/curvature metrics.
%  - Builds a feature table for all windows/signals.
%  - Runs PCA + kmeans clustering (k=2) to separate EMG-like vs EEG-like.
%  - Saves CSV, MAT and figures (spectrograms, mean spectra, PCA scatter, cluster summaries).
%
% Usage:
%  - Put this script in the same folder as your data.mat and run it in MATLAB.
%  - Check produced files: THz_features_table.csv, THz_results.mat, and PNG figures.
%
% NOTE: adjust fs (sampling frequency) if different from 256.

clc; clear; close all;
fprintf('THz_Style_Spectral_Analysis starting...\n');

%% --------- settings ----------
dataFile = 'data.mat';
if ~exist(dataFile,'file')
    error('data.mat not found in current folder.');
end
load(dataFile);    % load variables

% sampling frequency (change if needed)
fs = 256;

% analysis window for spectrogram
win_len_sec = 2;         % 2-second window
win_len = win_len_sec * fs;
noverlap = round(0.9 * win_len);  % 90% overlap (high time resolution)
nfft = 1024;

% pwelch params for mean spectrum estimation
pwelch_nperseg = 1024;

% signals to analyze (name, vector)
sig_list = {};
% add band signals if exist
if exist('dalt','var'), sig_list{end+1} = {'dalt', dalt(:)}; end
if exist('thalt','var'), sig_list{end+1} = {'thalt', thalt(:)}; end
if exist('alpha','var'), sig_list{end+1} = {'alpha', alpha(:)}; end
if exist('belta','var'), sig_list{end+1} = {'belta', belta(:)}; end

% raw signals
if exist('x1','var'), sig_list{end+1} = {'x1', x1(:)}; end
if exist('x2','var'), sig_list{end+1} = {'x2', x2(:)}; end

% if multichannel raw exists, add channel-average as representative
if exist('x12120x2Eedf','var')
    tmp = double(x12120x2Eedf);
    if ndims(tmp)==2
        % ensure samples x channels
        if size(tmp,1) < size(tmp,2), tmp = tmp'; end
        mean_sig = mean(tmp,2);
        sig_list{end+1} = {'x12120_mean', mean_sig(:)};
    end
end

nSignals = length(sig_list);
fprintf('Found %d signals to analyze.\n', nSignals);

%% create folders for outputs
outdir = 'THz_Spectral_Outputs';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% initialize storage
all_feature_rows = {}; % cell for table rows
feature_names = { ...
    'Signal','WindowStart','WindowEnd',...
    'totalPower','peakFreq','peakPower','specCentroid','specBandwidth','specFlatness','specEntropy',...
    'bp_delta','bp_theta','bp_alpha','bp_beta','bp_30_80','bp_80_200',...
    'RMS','MAV','ZCR','WL','kurtosis','spectrumCurvature','spectrumSlope' ...
    };

rowcount = 0;

%% main loop: per-signal compute spectrogram, mean spectrum and biomarkers
figure('visible','off'); % for safety when running on headless

for s = 1:nSignals
    name = sig_list{s}{1};
    sig = sig_list{s}{2};
    N = length(sig);
    fprintf('Analyzing signal %d/%d: %s (len=%d)\n', s, nSignals, name, N);
    
    % If signal is very long, we will process using sliding windows
    step = round(0.5 * win_len); % 50% step (overlap)
    starts = 1:step:(N - win_len + 1);
    if isempty(starts)
        warning('Signal %s too short for window length; using whole signal.', name);
        starts = 1;
        win_len_use = N;
    else
        win_len_use = win_len;
    end
    
    % create a figure for this signal spectrogram summary
    fig = figure('Name',['Spectrogram - ' name],'Visible','off','Units','normalized','Position',[0.1 0.1 0.8 0.6]);
    
    % We'll collect mean spectra across windows for FDA-like handling
    meanSpectra = [];
    freqs_for_pwelch = [];
    
    subplot_rows = ceil(length(starts)/2) + 1;
    spidx = 1;
    
    for k = 1:length(starts)
        istart = starts(k);
        iend = istart + win_len_use - 1;
        seg = sig(istart:iend);
        
        % Preprocessing: detrend and demean
        seg = detrend(seg);
        seg = seg - mean(seg);
        
        % Spectrogram (STFT)
        [S,F,T,P] = spectrogram(seg, hamming(win_len_use), round(0.9*win_len_use), nfft, fs);
        % S is complex STFT; P = |S|^2 power
        
        % Compute mean PSD from pwelch for this window (smoother estimate)
        % [Pxx,f_pxx] = pwelch(seg, pwelch_nperseg, round(0.5*pwelch_nperseg), nfft, fs);
        % --- safe pwelch parameters (ensure segment <= signal length) ---
        pwelch_nperseg_use = min(pwelch_nperseg, length(seg));          % don't ask pwelch for longer segment
        pwelch_noverlap = min(round(0.5*pwelch_nperseg_use), pwelch_nperseg_use-1);
        if pwelch_noverlap < 0, pwelch_noverlap = 0; end
        [Pxx,f_pxx] = pwelch(seg, pwelch_nperseg_use, pwelch_noverlap, nfft, fs);

        
        meanSpectra = [meanSpectra, Pxx]; %#ok<AGROW>
        freqs_for_pwelch = f_pxx;
        
        % Biomarkers from Pxx
        totalP = trapz(f_pxx, Pxx);
        [pkPow, idxpk] = max(Pxx);
        pkFreq = f_pxx(idxpk);
        centroid = sum(f_pxx .* Pxx) / (sum(Pxx) + eps);
        bandwidth = sqrt(sum(((f_pxx - centroid).^2) .* Pxx) / (sum(Pxx) + eps));
        % spectral flatness = geometric mean / arithmetic mean
        geoMean = exp(mean(log(Pxx + eps)));
        arithMean = mean(Pxx + eps);
        flatness = geoMean / arithMean;
        % spectral entropy (normalized)
        Pn = Pxx ./ (sum(Pxx)+eps);
        specEnt = -sum(Pn .* log2(Pn + eps)) / log2(length(Pn));
        
        % band powers
        bp_d = bandpower_from_psd(Pxx,f_pxx,[1 4]);
        bp_t = bandpower_from_psd(Pxx,f_pxx,[4 8]);
        bp_a = bandpower_from_psd(Pxx,f_pxx,[8 13]);
        bp_b = bandpower_from_psd(Pxx,f_pxx,[13 30]);
        bp_30_80 = bandpower_from_psd(Pxx,f_pxx,[30 80]);
        bp_80_200 = bandpower_from_psd(Pxx,f_pxx,[80 200]);
        
        % time-domain features
        RMSv = sqrt(mean(seg.^2));
        MAV = mean(abs(seg));
        ZCR = sum(abs(diff(sign(seg)))) / length(seg);
        WL = sum(abs(diff(seg)));
        K = kurtosis(seg);
        
        % FDA-like: smooth mean spectrum and compute slope/curvature
        % use smooth() on Pxx (moving average) - 'lowess' not necessary
        smoothP = smooth(f_pxx, Pxx, 0.05, 'loess');  % fraction 0.05, 'loess' available R2014a
        % derivative approx
        d1 = gradient(smoothP) ./ gradient(f_pxx+eps);
        d2 = gradient(d1) ./ gradient(f_pxx+eps);
        specCurv = sum(abs(d2));
        % spectrum slope = linear regression of log(Pxx) vs log(f)
        % consider frequency range > 1 Hz to avoid log(0)
        valid_idx = f_pxx > 1;
        if sum(valid_idx) > 2
            p = polyfit(log(f_pxx(valid_idx)), log(Pxx(valid_idx)+eps), 1);
            slope = p(1);
        else
            slope = NaN;
        end
        
        % store a row
        rowcount = rowcount + 1;
        all_feature_rows{rowcount,1} = name;
        all_feature_rows{rowcount,2} = istart;
        all_feature_rows{rowcount,3} = iend;
        all_feature_rows{rowcount,4} = totalP;
        all_feature_rows{rowcount,5} = pkFreq;
        all_feature_rows{rowcount,6} = pkPow;
        all_feature_rows{rowcount,7} = centroid;
        all_feature_rows{rowcount,8} = bandwidth;
        all_feature_rows{rowcount,9} = flatness;
        all_feature_rows{rowcount,10} = specEnt;
        all_feature_rows{rowcount,11} = bp_d;
        all_feature_rows{rowcount,12} = bp_t;
        all_feature_rows{rowcount,13} = bp_a;
        all_feature_rows{rowcount,14} = bp_b;
        all_feature_rows{rowcount,15} = bp_30_80;
        all_feature_rows{rowcount,16} = bp_80_200;
        all_feature_rows{rowcount,17} = RMSv;
        all_feature_rows{rowcount,18} = MAV;
        all_feature_rows{rowcount,19} = ZCR;
        all_feature_rows{rowcount,20} = WL;
        all_feature_rows{rowcount,21} = K;
        all_feature_rows{rowcount,22} = specCurv;
        all_feature_rows{rowcount,23} = slope;
        
        % plot spectrogram of first few windows into figure
        if k <= 4
            subplot(3,4,spidx);
            imagesc(T, F, 10*log10(abs(P))); axis xy;
            title(sprintf('%s window %d (t:%d-%d s)', name, k, round(istart/fs), round(iend/fs)));
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            colormap jet; colorbar;
            spidx = spidx + 1;
        end
    end % windows loop
    
    % plot mean spectrum (average of Pxx across windows) in a separate axis
    meanP = mean(meanSpectra,2);
    subplot(3,4,spidx);
    plot(freqs_for_pwelch, meanP, 'LineWidth', 1.4);
    xlabel('Freq (Hz)'); ylabel('Power');
    title([name ' - mean spectrum (averaged windows)']);
    grid on;
    
    % save figure per signal
    saveas(fig, fullfile(outdir, sprintf('Spectrogram_and_mean_%s.png', name)));
    close(fig);
end % signals

%% Build table and save CSV/MAT
% convert all_feature_rows to table
nRows = size(all_feature_rows,1);
T = cell2table(all_feature_rows, 'VariableNames', feature_names);
writetable(T, fullfile(outdir,'THz_features_table.csv'));
save(fullfile(outdir,'THz_results.mat'), 'T', 'feature_names');

fprintf('Feature extraction done. Saved features CSV and MAT in folder %s\n', outdir);

%% PCA and clustering on features
% choose numeric columns 4:end
Xfeat = table2array(T(:,4:end));
% replace NaN with 0 for PCA (or mean impute)
Xfeat(isnan(Xfeat)) = 0;
% zscore
muX = mean(Xfeat,1); sigmaX = std(Xfeat,0,1) + eps;
Xz = (Xfeat - repmat(muX,size(Xfeat,1),1)) ./ repmat(sigmaX,size(Xfeat,1),1);

% PCA
[coeff, score, latent, tsquared, explained] = pca(Xz);

% kmeans clustering (k=2)
k = 2;
rng(2);
[idx_km, Ckm] = kmeans(Xz, k, 'Replicates', 20);

% decide which cluster is EMG-like by checking mean bp_30_80 (feature index 15 in table)
col_bp30 = 15 - 3; % because Xfeat starts at col 4 in T -> index 11 in Xfeat? careful: direct mapping easier:
% compute cluster means on original T field bp_30_80
bp30_all = T.bp_30_80;
clusterMeans = zeros(k,1);
for c=1:k
    clusterMeans(c) = mean(bp30_all(idx_km==c));
end
[~, emg_cluster] = max(clusterMeans);
eeg_cluster = setdiff(1:k, emg_cluster);

% add clustering results into table and save
T.Cluster = idx_km;
T.Is_EMG_like = (idx_km == emg_cluster);

writetable(T, fullfile(outdir,'THz_features_table_with_clusters.csv'));
save(fullfile(outdir,'THz_results_with_clusters.mat'), 'T', 'coeff', 'score', 'explained', 'idx_km', 'emg_cluster', 'eeg_cluster');

%% plot PCA scatter colored by cluster
fig = figure('Visible','off');
gscatter(score(:,1), score(:,2), idx_km);
xlabel('PC1'); ylabel('PC2'); title('PCA scatter colored by kmeans cluster');
legend('Cluster 1','Cluster 2');
saveas(fig, fullfile(outdir,'THz_PCA_clusters.png'));
close(fig);

%% Summarize cluster assignments per signal
uniqueSignals = unique(T.Signal);
summary = cell(length(uniqueSignals),3);
for i=1:length(uniqueSignals)
    nm = uniqueSignals{i};
    mask = strcmp(T.Signal, nm);
    total = sum(mask);
    emgw = sum(T.Is_EMG_like(mask));
    summary{i,1} = nm;
    summary{i,2} = total;
    summary{i,3} = emgw;
    fprintf('%s: windows=%d, EMG_like_windows=%d (%.1f%%)\n', nm, total, emgw, 100*emgw/total);
end

% save summary as CSV
summaryTable = cell2table(summary, 'VariableNames', {'Signal','TotalWindows','EMG_like_Windows'});
writetable(summaryTable, fullfile(outdir,'THz_signal_cluster_summary.csv'));
save(fullfile(outdir,'THz_signal_cluster_summary.mat'), 'summaryTable');

fprintf('Clustering done. Summary saved. Check folder: %s\n', outdir);
fprintf('THz-style spectral analysis finished.\n');

%% -------------- helper functions --------------
% function bp = bandpower_from_psd(Pxx,f,band)
% function H = spectral_entropy(Pxx)

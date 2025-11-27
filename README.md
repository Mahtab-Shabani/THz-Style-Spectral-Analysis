# THz-Style-Spectral-Analysis
THz-style Spectral Biomarker Extraction for EEG/EMG signals
This MATLAB script performs a time-frequency and spectral biomarker analysis inspired by THz spectroscopy approaches. It is compatible with MATLAB R2014a.

<img width="1200" height="900" alt="THz_PCA_clusters" src="https://github.com/user-attachments/assets/62187465-9a8d-475b-9825-b51ff8d41424" />

## Features
* Loads data.mat (expects: dalt, thalt, alpha, belta, x1, x2, x12120x2Eedf)
* Computes per-signal spectrograms (STFT) and mean spectra (via pwelch)
* Extracts spectral biomarkers:
    1. Entropy, centroid, bandwidth, flatness, peak frequency & power
    2. Band powers (delta, theta, alpha, beta, 30–80 Hz, 80–200 Hz)
    3. Time-domain features (RMS, MAV, ZCR, waveform length, kurtosis)
    4. FDA-style metrics (smoothed spectrum, slope, curvature)
* Builds a feature table for all signals and windows
* Performs PCA + k-means clustering (k=2) to separate EMG-like vs EEG-like signals
* Saves results as CSV, MAT files, and figures (spectrograms, mean spectra, PCA scatter, cluster summaries)

## Usage
1. Place THz_Style_Spectral_Analysis.m in the same folder as your data.mat
2. Open MATLAB and run the script
3. Outputs will be saved in the THz_Spectral_Outputs folder:
    * THz_features_table.csv
    * THz_results.mat
    * THz_features_table_with_clusters.csv
    * PNG figures (Spectrogram, Mean Spectrum, PCA scatter)

## Notes
* Default sampling frequency: fs = 256 Hz (change if needed)
* Window length: 2 s with 90% overlap for high temporal resolution
* Works with single-channel or multichannel signals (averaged if multichannel)

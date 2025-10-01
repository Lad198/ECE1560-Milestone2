%% 1560 Milestone 2 - Logan Duncan - Time-Domain %%
clear; close all; clc;
%% ---------------------- USER INPUT (choose ONE) ----------------------
% A) File is in the SAME folder as this script:
%wavFile = 'cardinal_chirp_single.wav';
% B) Or, hardcode the FULL path:
wavFile = 'C:\Users\logan\Downloads\ECE1560_fall2025_project_ms2_audio\alto_sax\alto_sax_a4.wav';
% Validate the path before doing anything else
if ~exist(wavFile,'file')
    error('File not found: %s', wavFile);
end
% Read the audio file
[x, fs] = audioread(wavFile);
% ---- Guard against empty/degenerate files ----
if isempty(x)
    error('Audio file contains zero samples: %s', wavFile);
end
% ---- Handle multi-channel audio (robust for unknown files) ----
% For the possibility sterio-audio is given as input this function
% Converts to mono by averaging channels to ensure a single trace for the plot.
if size(x,2) > 1
    x = mean(x, 2);
    channelInfo = ' (mono mixdown)';
else
    channelInfo = '';
end
% ---- Basic parameters ----
N  = length(x);   % number of samples
Ts = 1/fs;        % time between samples [s]
t  = (0:N-1) * Ts; % time vector from 0 to (N-1)*Ts
%% ---------------------- ACF EXPLORATION (not the final algorithm) ----------------------
F0Range   = [70 500];   % Hz: tune to your instrument set
MinProm   = 0.04;       % Peak prominence for display (helps suppress tiny bumps)
DoClip    = false;      % Try true later to see how center-clipping changes ACF
ClipAlpha = 0.6;        % 0..1 aggressiveness if DoClip = true
% 1) Light preprocess (keeps ACF cleaner)
x = x - mean(x);                             % remove DC
w = ones(N,1);              % rectangular window (no taper)
xw = x .* w;                % unchanged
if DoClip
    % Optional experiment: center clipping emphasizes periodic fine structure
    Tclip = ClipAlpha * max(abs(x));
    xc = zeros(size(x));
    pos = x >  Tclip;  neg = x < -Tclip;
    xc(pos) = x(pos) - Tclip;
    xc(neg) = x(neg) + Tclip;
    xw = xc .* w;                             % analyze the clipped version
end
% 2) Normalized ACF for nonnegative lags (simple xcorr path)
[r_full, lags_full] = xcorr(xw);             % length 2N-1
keep = lags_full >= 0;                       % keep k >= 0
r = r_full(keep);
lags = lags_full(keep);                      % in samples (0..N-1)
% Normalize so R(0)=1 (makes peaks comparable across files)
if r(1) > 0
    r = r / r(1);
else
    % very rare fallback
    r = r - min(r); r = r / max(r + eps);
end
% 3) Map f0 range to a lag (sample) range and inspect peaks
kmin = max(2, floor(fs / F0Range(2)));       % skip lag=0 (trivial peak)
kmax = min(numel(r)-1, ceil(fs / F0Range(1)));
seg = r(kmin:kmax);
tau = (kmin:kmax) / fs;                      % seconds for those lags
% Find local peaks inside the candidate region (for inspection only)
isPk = islocalmax(seg, 'MinProminence', MinProm);
pk_idx_local = find(isPk);
pk_tau  = tau(isPk);                         % seconds to first few peaks
pk_vals = seg(isPk);
% Rough observed spacing between successive peaks (useful sanity check)
if numel(pk_tau) >= 3
    dT = diff(pk_tau);
    T0_obs = median(dT);
    f0_obs = 1 / T0_obs;
else
    T0_obs = NaN; f0_obs = NaN;
end
% ---- Print observations you can cite in your write-up ----
fprintf('\n--- ACF exploration for: %s%s ---\n', wavFile, channelInfo);
fprintf('Fs = %d Hz, N = %d (%.2f s). ACF normalized: R(0)=1.\n', fs, N, N/fs);
fprintf('Search window: k=[%d..%d] -> tau=[%.5f..%.5f] s (f0≈[%d..%d] Hz)\n', ...
        kmin, kmax, kmin/fs, kmax/fs, F0Range(1), F0Range(2));
fprintf('Found %d local peaks in window.\n', numel(pk_tau));
if ~isnan(f0_obs)
    fprintf('Observed inter-peak spacing T0≈%.5f s -> f0≈%.2f Hz (exploratory)\n', T0_obs, f0_obs);
end
for i = 1:min(5, numel(pk_tau))
    fprintf('  Peak %d: tau=%.5f s (~%.1f Hz), value=%.3f\n', i, pk_tau(i), 1/pk_tau(i), pk_vals(i));
end
% ---- Plots: raw snippet + normalized ACF with peaks and shaded search region ----
ShowSecondsOfWave = 0.12; ns = min(N, round(ShowSecondsOfWave*fs));
figure('Name','Waveform + ACF'); 
subplot(2,1,1);
plot(t(1:ns), x(1:ns)); grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Waveform snippet (%d ms)%s', round(1000*ns/fs), channelInfo));
subplot(2,1,2);
lag_s = lags / fs;
plot(lag_s, r, 'LineWidth', 1); hold on; grid on;
yl = ylim;
patch([kmin kmax kmax kmin]/fs, [yl(1) yl(1) yl(2) yl(2)], ...
      [0.9 0.95 1], 'EdgeColor','none', 'FaceAlpha', 0.2);
if ~isempty(pk_tau)
    stem(pk_tau, r(kmin + pk_idx_local - 1), 'filled', 'MarkerSize', 4);
end
xlabel('Lag \tau (s)'); ylabel('R_x(\tau) / R_x(0)');
ttl = 'Normalized ACF (nonnegative lags)';
if DoClip, ttl = [ttl, ' — center-clipped']; end
title(ttl);
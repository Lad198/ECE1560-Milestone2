filename = "alto sax fsharp4.wav"
[y, Fs] = audioread(filename);
% Use the entire signal
x = y;
N = length(x);
% Remove DC offset
x = x - mean(x);
% Compute FFT
X = fft(x);
mag = abs(X(1:floor(N/2)));  % Single-sided magnitude spectrum
f = (0:floor(N/2)-1) * Fs / N;  % Frequency vector
% Find peaks in the magnitude spectrum
[pks, locs] = findpeaks(mag, 'MinPeakHeight', max(mag) * 0.1);
% Handle case with no peaks
if isempty(locs)
f0 = 0;
fprintf('No significant peaks found.\n');
return;
end
% Assume the first significant peak is the fundamental
f0 = f(locs(1))
%% 1560 Milestone 2 - TD_A Batch Compare (Normalized ACF vs Excel)
clear; close all; clc;

%% ---------------------- USER INPUT ----------------------
baseDir   = 'C:\Users\logan\Downloads\ECE1560_fall2025_project_ms2_audio (1)'; % parent folder with instrument subfolders
excelPath = '';   % if you have the Excel/CSV of expected values, put path here ('.xlsx' or '.csv'); leave empty to use embedded table
PctTol    = 5;    % ±5% threshold
F0Range   = [70 500]; % Hz
MinProm   = 0.03; % ACF peak prominence

% filename pitch tokens to detect (lowercase)
validPitches = ["a3","bflat3","b3","c4","csharp4","d4","eflat4","e4","f4","fsharp4","g4","gsharp4","a4"];

%% ---------------------- Load reference (Excel or embedded) ----------------------
if ~isempty(excelPath) && isfile(excelPath)
    [~,~,ext] = fileparts(excelPath);
    switch lower(ext)
        case '.xlsx'
            T = readtable(excelPath);
        case '.csv'
            T = readtable(excelPath, 'TextType','string');
        otherwise
            error('Unsupported reference file: %s', excelPath);
    end
    % Try to normalize column names if they differ
    T.Properties.VariableNames = lower(strrep(T.Properties.VariableNames,' ','_'));
    % Expect columns: instrument, pitch, f0_ideal (or f0)
    if ismember('f0', T.Properties.VariableNames) && ~ismember('f0_ideal', T.Properties.VariableNames)
        T.f0_ideal = T.f0; 
    end
    refT = table(string(lower(string(T.instrument))), string(lower(string(T.pitch))), T.f0_ideal, ...
                 'VariableNames',{'instrument','pitch','f0_ideal'});
else
    refT = buildEmbeddedReference(); % from your pasted values
end

%% ---------------------- Walk files ----------------------
files = dir(fullfile(baseDir,'**','*.wav'));
if isempty(files), error('No WAV files found under: %s', baseDir); end
fprintf('Found %d wav files.\n', numel(files));

Results = table();
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    [parentDir, fnameNoExt, ~] = fileparts(fpath);
    [~, instrument] = fileparts(parentDir);
    instrument = string(lower(instrument));
    pitch = extractPitchToken(lower(string(fnameNoExt)), validPitches);
    if pitch==""
        warning('No pitch token in filename, skipping: %s', fpath);
        continue;
    end

    % lookup ideal
    row = refT(refT.instrument==instrument & refT.pitch==pitch, :);
    if height(row)~=1
        warning('No reference for %s / %s, skipping: %s', instrument, pitch, fpath);
        continue;
    end
    f0_ideal = row.f0_ideal;

    % --- Read audio & TD_A estimate ---
    [x, fs] = audioread(fpath);
    if isempty(x), warning('Empty audio: %s', fpath); continue; end
    if size(x,2) > 1, x = mean(x,2); end
    x = x - mean(x);

    % take loudest ~0.4 s (avoid silence/attack bias)
    x = takeLoudestWindow(x, fs, 0.40, 0.10);

    % Normalized ACF, nonzero-lag strongest peak with parabolic interpolation
    f0_est = f0_from_acf_peak(x, fs, F0Range, MinProm);

    % Metrics
    cents_err = 1200 * log2(f0_est / f0_ideal);
    pct_err   = 100  * (f0_est - f0_ideal) / f0_ideal;
    pass      = abs(pct_err) <= PctTol;

    Results = [Results; table(string(fpath), instrument, pitch, fs, f0_ideal, f0_est, cents_err, pct_err, pass)]; %#ok<AGROW>
end

Results.Properties.VariableNames = {'file','instrument','pitch','fs','f0_ideal','TD_A_f0','TD_A_cents','TD_A_pct','TD_A_pass'};
writetable(Results, 'results_td_a.csv');
fprintf('Saved TD_A results to results_td_a.csv\n');

% Quick summary
acc = 100*mean(Results.TD_A_pass);
medC = median(abs(Results.TD_A_cents),'omitnan');
fprintf('[TD_A] Accuracy: %.1f%% pass @ ±%d%%, median |cents| = %.1f\n', acc, PctTol, medC);

G = Results(:,{'instrument','TD_A_cents','TD_A_pass'});
G.absC = abs(G.TD_A_cents);
S = groupsummary(G,'instrument',{'mean','median'},{'TD_A_pass','absC'});
S.mean_TD_A_pass = 100*S.mean_TD_A_pass;
S.Properties.VariableNames{'median_absC'} = 'median_abs_cents';
disp(sortrows(S(:,{'instrument','mean_TD_A_pass','median_abs_cents'}),'median_abs_cents'));

%% ============================ Helpers ============================
function f0 = f0_from_acf_peak(x, fs, F0Range, MinProm)
    [xc,lags] = xcorr(x,'coeff');
    keep = lags>=0;
    r = xc(keep);
    % map F0 to lag range
    kmin = max(2, floor(fs / F0Range(2)));
    kmax = min(numel(r)-1, ceil(fs / F0Range(1)));
    if kmax <= kmin+1, f0 = NaN; return; end
    seg = r(kmin:kmax);
    [pkVals, pkLocs] = findpeaks(seg, 'MinPeakProminence', MinProm);
    if isempty(pkLocs), f0 = NaN; return; end
    [~,idx] = max(pkVals);
    L = kmin + pkLocs(idx) - 1;

    % parabolic interpolation about L
    if L>1 && L<numel(r)
        y1=r(L-1); y2=r(L); y3=r(L+1);
        denom = y1 - 2*y2 + y3;
        delta = 0;
        if denom~=0
            delta = 0.5*(y1 - y3)/denom; % [-0.5,0.5] samples
        end
    else
        delta = 0;
    end
    Lhat = L + delta;
    f0 = fs / Lhat;
end

function xw = takeLoudestWindow(x, fs, segS, hopS)
    win = max(1, round(segS*fs));
    if numel(x) <= win, xw = x; return; end
    hop = max(1, round(hopS*fs));
    idxs = 1:hop:(numel(x)-win+1);
    rmsVals = zeros(numel(idxs),1);
    for i = 1:numel(idxs)
        ii = idxs(i);
        rmsVals(i) = sqrt(mean(x(ii:ii+win-1).^2));
    end
    [~,m] = max(rmsVals);
    ii = idxs(m);
    xw = x(ii:ii+win-1);
end

function p = extractPitchToken(fname, validPitches)
    p = "";
    for tok = validPitches
        if contains(fname, tok)
            p = tok; return;
        end
    end
end

function refT = buildEmbeddedReference()
% (instrument, pitch, f0_ideal) from the instructor list you pasted
rows = {
"alto_sax","a3",220
"alto_sax","bflat3",233.0818808
"alto_sax","b3",246.9416506
"alto_sax","c4",261.6255653
"alto_sax","csharp4",277.182631
"alto_sax","d4",293.6647679
"alto_sax","eflat4",311.1269837
"alto_sax","e4",329.6275569
"alto_sax","f4",349.2282314
"alto_sax","fsharp4",369.9944227
"alto_sax","g4",391.995436
"alto_sax","gsharp4",415.3046976
"alto_sax","a4",440
"clarinet","a3",220
"clarinet","bflat3",233.0818808
"clarinet","b3",246.9416506
"clarinet","c4",261.6255653
"clarinet","csharp4",277.182631
"clarinet","d4",293.6647679
"clarinet","eflat4",311.1269837
"clarinet","e4",329.6275569
"clarinet","f4",349.2282314
"clarinet","fsharp4",369.9944227
"clarinet","g4",391.995436
"clarinet","gsharp4",415.3046976
"clarinet","a4",440
"piano","a3",220
"piano","bflat3",233.0818808
"piano","b3",246.9416506
"piano","c4",261.6255653
"piano","csharp4",277.182631
"piano","d4",293.6647679
"piano","eflat4",311.1269837
"piano","e4",329.6275569
"piano","f4",349.2282314
"piano","fsharp4",369.9944227
"piano","g4",391.995436
"piano","gsharp4",415.3046976
"piano","a4",440
"tenor_sax","a3",220
"tenor_sax","bflat3",233.0818808
"tenor_sax","b3",246.9416506
"tenor_sax","c4",261.6255653
"tenor_sax","csharp4",277.182631
"tenor_sax","d4",293.6647679
"tenor_sax","eflat4",311.1269837
"tenor_sax","e4",329.6275569
"tenor_sax","f4",349.2282314
"tenor_sax","fsharp4",369.9944227
"tenor_sax","g4",391.995436
"tenor_sax","gsharp4",415.3046976
"tenor_sax","a4",440
"trumpet","a3",220
"trumpet","bflat3",233.0818808
"trumpet","b3",246.9416506
"trumpet","c4",261.6255653
"trumpet","csharp4",277.182631
"trumpet","d4",293.6647679
"trumpet","eflat4",311.1269837
"trumpet","e4",329.6275569
"trumpet","f4",349.2282314
"trumpet","fsharp4",369.9944227
"trumpet","g4",391.995436
"trumpet","gsharp4",415.3046976
"trumpet","a4",440
"violin","a3",220
"violin","bflat3",233.0818808
"violin","b3",246.9416506
"violin","c4",261.6255653
"violin","csharp4",277.182631
"violin","d4",293.6647679
"violin","eflat4",311.1269837
"violin","e4",329.6275569
"violin","f4",349.2282314
"violin","fsharp4",369.9944227
"violin","g4",391.995436
"violin","gsharp4",415.3046976
"violin","a4",440
"voice","a3",220
"voice","bflat3",233.0818808
"voice","b3",246.9416506
"voice","c4",261.6255653
"voice","csharp4",277.182631
"voice","d4",293.6647679
"voice","eflat4",311.1269837
"voice","e4",329.6275569
"voice","f4",349.2282314
"voice","fsharp4",369.9944227
"voice","g4",391.995436
"voice","gsharp4",415.3046976
"voice","a4",440
};
refT = cell2table(rows, 'VariableNames',{'instrument','pitch','f0_ideal'});
refT.instrument = string(lower(refT.instrument));
refT.pitch      = string(lower(refT.pitch));
end

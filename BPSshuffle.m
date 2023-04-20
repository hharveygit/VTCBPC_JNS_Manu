%% Demonstrate that BPSs are invariant to temporal/spatial order
%
% 2023/04
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Creates spectrograms and broadband traces for consensus BPCs
%    Copyright (C) 2022  Harvey Huang
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% Configure user parameters, paths

subNames = {'3'}; % sub3 presented in
subChs = {{'LB6'}};
hemis = 'l';

%% Load and concatenate voltage BPC outputs and spectrograms

specData = cell(0, 4);
srate = 2048;

sub = subNames{1};
ch = subChs{1}{1};
    
fprintf('sub-%s, ch-%s\n', sub, ch);

% Load spectrogram data
load(fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_spectralBroadband.mat', sub, ch)));

S = log10(S(f >= 12, :, :)); % keep only frequencies higher than 12 Hz to reduce dimensionality. Convert to log power
fspec = f(f >= 12);

Sdown = zeros(length(fspec), size(S, 2)/16, size(S, 3));
for kk = 1:length(fspec)
    Sdown(kk, :, :) = resample(squeeze(double(S(kk, :, :))), 1, 16); % downsample by factor of 16
end
ttDown = linspace(tt(1), tt(end), size(S, 2)/16); % downsampled time vector
Sdown = Sdown - mean(Sdown(:, ttDown > -0.5 & ttDown < -0.1, :), 2); % log-normalize spectrograms by baseline, far from stim
Sdown = Sdown(fspec >= 12, ttDown >= 0.1 & ttDown <= 0.5, :);

Sdownvec = reshape(Sdown, size(Sdown, 1)*size(Sdown, 2), size(Sdown, 3)); % flatten first dimension
specData = [specData; {sub}, {ch}, {events_ch}, {Sdownvec}];

% Load voltage BPC stuff
snrTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCSNR.tsv', sub, ch)), 'electrical_stimulation_site', 1);  
stimSites = snrTable.electrical_stimulation_site;

clear S f sigdataBB tt pxx pxxWinds fPSD

%% Plot the spectrograms to go into the BPS algorithm

ttBPC = ttDown(ttDown >= 0.1 & ttDown <= 0.5);
cmapSpectra = getCmapSpec();

% Sdown is F by time by trial
for tr = 1:30:height(events_ch)
    figure; uimagesc(ttBPC, fspec, Sdown(:, :, tr), [-2, 2]);
    colormap(cmapSpectra);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    axis xy;
    colorbar;
    saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_spec%02d', specData{1, 1}, specData{1, 2}, tr)), 'svg');
    close(gcf)
end

% plot the flattened form
figure('Position', [200, 200, 400, 1000]);
imagesc(Sdownvec, [-2, 2]);
colormap(cmapSpectra);
colorbar;
saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_flattenedSpecs', specData{1, 1}, specData{1, 2})), 'svg');
close(gcf)


%% Calculate basis profile spectrograms

mkdir('output/BPS/shuffled');

rng('default'); % for reproducibility

ttBPC = ttDown(ttDown >= 0.1 & ttDown <= 0.5);

ii = 1;

events_ch = specData{ii, 3};
sites = groupby(events_ch.electrical_stimulation_site);
sites(cellfun(@length, sites(:, 2)) < 8, :) = []; % exclude sites with < 8 trials
pairTypes = struct('pair', sites(:, 1), 'indices', sites(:, 2));

[B, exc, ~] = bpc_identify(Sdownvec, pairTypes, 1, 50); % 50 reruns, log10 psd
BPSflat = [B.curve];

% plot the flat BPSs
f = figure('Position',[100 100 400 600]); hold on
set(gca, 'YTicklabel', []);
plotTrials(1:length(BPSflat), BPSflat, 0.1, [], [], 'LineWidth', 1);
saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_specBPCFlat%d', specData{ii, 1}, specData{ii, 2})), 'svg');

% pull out spectrograms and reshape
Bspecs = zeros(length(fspec), length(ttBPC), length(B));
for bb = 1:length(B)
    Bspecs(:, :, bb) = reshape(B(bb).curve, length(fspec), length(ttBPC));
end

% sort BPCs such that most negative BPC is first, for plotting consistency across electrodes & subjects
[~, bpcOrd] = sort(squeeze(sum(Bspecs(fspec >= 60 & fspec <= 120, ttBPC >= 0.15 & ttBPC <= 0.35, :), [1, 2])));
Bspecs = Bspecs(:, :, bpcOrd);
B = B(bpcOrd);

% reshape each spectrogram and plot
for bb = 1:length(B)

    figure; uimagesc(ttBPC, fspec, Bspecs(:, :, bb), [-0.03, 0.03]);
    colormap(cmapSpectra);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    axis xy;
    colorbar;
    saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_specBPC%d', specData{ii, 1}, specData{ii, 2}, bb)), 'png');
    saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_specBPC%d', specData{ii, 1}, specData{ii, 2}, bb)), 'svg');
    close(gcf)

end

%% Shuffle spectrograms first and then calculate BPS

% reshuffle the timebins and plot
rng('default'); % for reproducibility
order = randperm(length(Sdownvec));
SdownvecShuffled = Sdownvec(order, :);
figure('Position', [200, 200, 400, 1000]);
imagesc(SdownvecShuffled, [-2, 2]);
colormap(cmapSpectra);
colorbar;
saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_flattenedSpecsShuffled', specData{1, 1}, specData{1, 2})), 'svg');
close(gcf);

[B, exc, ~] = bpc_identify(SdownvecShuffled, pairTypes, 1, 50); % 50 reruns, log10 psd

% plot the flat BPSs (shuffled)
BPSflatShuffled = [B.curve];
f = figure('Position',[100 100 400 600]); hold on
set(gca, 'YTicklabel', []);
plotTrials(1:length(BPSflatShuffled), BPSflatShuffled, 0.1, [], [], 'LineWidth', 1);
saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_specBPCFlatShuffled%d', specData{ii, 1}, specData{ii, 2})), 'svg');

% unshuffle flat BPSs and plot
[~, idx] = sort(order); % sort the permutation order
BPSflatUnshuffled = BPSflatShuffled(idx, :);
f = figure('Position',[100 100 400 600]); hold on
set(gca, 'YTicklabel', []);
plotTrials(1:length(BPSflatUnshuffled), BPSflatUnshuffled, 0.1, [], [], 'LineWidth', 1);
saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_specBPCFlatUnshuffled%d', specData{ii, 1}, specData{ii, 2})), 'svg');

% pull out spectrograms and reshape
Bspecs = zeros(length(fspec), length(ttBPC), length(B));
for bb = 1:length(B)
    Bspecs(:, :, bb) = reshape(BPSflatUnshuffled(:, bb), length(fspec), length(ttBPC));
end

% sort BPCs such that most negative BPC is first, for plotting consistency across electrodes & subjects
[~, bpcOrd] = sort(squeeze(sum(Bspecs(fspec >= 60 & fspec <= 120, ttBPC >= 0.15 & ttBPC <= 0.35, :), [1, 2])));
Bspecs = Bspecs(:, :, bpcOrd);
B = B(bpcOrd);

% reshape each spectrogram and plot
for bb = 1:length(B)
  
    figure; uimagesc(ttBPC, fspec, Bspecs(:, :, bb), [-0.03, 0.03]);
    colormap(cmapSpectra);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    axis xy;
    colorbar;
    saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_specBPCUnshuffled%d', specData{ii, 1}, specData{ii, 2}, bb)), 'png');
    saveas(gcf, fullfile('output', 'BPS', 'shuffled', sprintf('sub-%s_ch-%s_specBPCUnshuffled%d', specData{ii, 1}, specData{ii, 2}, bb)), 'svg');
    close(gcf)

end

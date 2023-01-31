%% File for calculating spectra and broadband for each consensus BPC
% 2021/10/20
% use after clusterBPC determines consensus BPCs in BPCSnr
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
%% Configure subjects. Load PSD, spectrograms, and broadband signals. Normalize by baseline, store trial means into cell array

clear;
set(0, 'DefaultFigureRenderer', 'painters');
cmapSpectra = getCmapSpec();

subNames = {'1', '2', '3', '4'}; % Collateral sulcus measurement electrodes only
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6'}};

masterData = cell(0, 7); % sub, ch, consBPC, SNR, stimsite, spectrogram, sigdataBB
for ii = 1:length(subNames)
    sub = subNames{ii};
    for jj = 1:length(subChs{ii})
        ch = subChs{ii}{jj};
        
        load(fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_spectralBroadband.mat', sub, ch))); % S, f, tt, sigdataBB
        snrTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCSNR.tsv', sub, ch)), 'electrical_stimulation_site', 1);
        
        stimSites = snrTable.electrical_stimulation_site;
        for kk = 1:length(stimSites)
            
            Scurr = S(:, :, strcmp(events_ch.electrical_stimulation_site, stimSites{kk}));
            Scurr = Scurr./geomean(Scurr(:, tt > -0.5 & tt < -0.1, :), 2); % normalize by baseline, staying far away from stim artifact
            Scurr = geomean(Scurr, 3); % average across trials
            
            BBcurr = sigdataBB(:, strcmp(events_ch.electrical_stimulation_site, stimSites{kk}));
            BBcurr = BBcurr.^2; % convert to power
            BBcurr = BBcurr./geomean(BBcurr(tt > -0.5 & tt <- 0.1, :)); % normalize by baseline
            BBcurr = geomean(BBcurr, 2); % mean across trials
            
            fprintf('sub-%s, ch-%s, stim-%s\n', sub, ch, stimSites{kk});
            masterData = [masterData; {sub}, {ch}, {snrTable.consBPC(kk)}, {snrTable.SNR(kk)}, stimSites(kk), {Scurr}, {BBcurr}];
        end
        
    end
end

%% Save example spectrograms from individual subjects for first consensus BPC (for figure 5A, bottom)

mkdir(fullfile('output', 'spectralBB', 'spectrogramsConsBPC1'));

for ii = 1:size(masterData, 1)
    if ~(masterData{ii, 3} == 1) || masterData{ii, 4} < 1, continue; end % only save examples for first consensus BPC sites with significant SNR
    
    fig = figure;
    uimagesc(tt, f, log10(masterData{ii, 6}), [-3, 3]);
    colormap(cmapSpectra);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xlim([-0.1, 1]);
    axis xy;
    colorbar;
    
    saveas(fig, fullfile('output', 'spectralBB', 'spectrogramsConsBPC1', sprintf('spectrogram_sub-%s_ch-%s_stim-%s', masterData{ii, 1}, masterData{ii, 2}, masterData{ii, 5})), 'svg');
    close(fig);
end

%% Calculate t-statistic spectrogram heatmaps for each consensus BPC. Visualize

tstatSpecs = zeros(length(f), length(tt), 5); % hold tstat for each consensus BPC. #5 = excluded CCEPs
pSpecs = zeros(size(tstatSpecs)); % hold p-value for each consensus BPC

figure('Position', [600, 200, 600, 1600]);
for ii = 1:4
    fprintf('Consensus BPC %d\n', ii);
    specsCurr = cat(3, masterData{[masterData{:, 3}] == ii & [masterData{:, 4}] >= 1, 6}); % spectrograms for current BPC w SNR > 1

    for jj = 1:length(f)
        fprintf('%.1fHz, ', f(jj));
        for kk = 1:length(tt)
            [~, p, ~, stats] = ttest(log10(squeeze(specsCurr(jj, kk, :)))); % ttest at each time-freq bin
            tstatSpecs(jj, kk, ii) = stats.tstat;
            pSpecs(jj, kk, ii) = p;
        end
    end
    fprintf('\n');

    subplot(5, 1, ii);
    uimagesc(tt, f, tstatSpecs(:, :, ii), [-10, 10]);
    colormap(cmapSpectra);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xlim([-0.1, 1]);
    ylabel('t-statistic');
    axis xy;
    colorbar;
    title(sprintf('consensus BPC %d', ii));
end

disp('Excluded CCEPs');
specsExc = cat(3, masterData{isnan([masterData{:, 3}]), 6});
for jj = 1:length(f)
    fprintf('%.1fHz, ', f(jj));
    for kk = 1:length(tt)
        [~, p, ~, stats] = ttest(log10(squeeze(specsExc(jj, kk, :))));
        tstatSpecs(jj, kk, 5) = stats.tstat; % excluded cceps are treated as "5th" consensus BPC
        pSpecs(jj, kk, 5) = p;
    end
end
fprintf('\n');

subplot(5, 1, 5);
uimagesc(tt, f, tstatSpecs(:, :, 5), [-10, 10]);
colormap(cmapSpectra);
xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xlim([-0.1, 1]);
ylabel('t-statistic');
axis xy;
colorbar;
title('Excluded CCEPs');
xlabel('time (s)');

%% Save consensus spectrograms individually to file

xlims = [-0.1, 1];
ttLim = tt(tt>=xlims(1) & tt<xlims(2));

% consensus BPC "5" is the group of stim sites excluded by BPC algorithm
for ii = 1:5
    fig = figure;
    uimagesc(ttLim, f, tstatSpecs(:, tt>=xlims(1) & tt <=xlims(2), ii), [-10, 10]);
    colormap(cmapSpectra);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ylabel('t-statistic');
    axis xy;
    colorbar;
    if ii == 5
        xlabel('time (s)');
        title('Excluded CCEPs');
    else
        title(sprintf('consensus BPC %d', ii));
    end

    saveas(fig, fullfile('output', 'spectralBB', sprintf('spectrogram_consBPCs_num%d', ii)), 'svg');
    saveas(fig, fullfile('output', 'spectralBB', sprintf('spectrogram_consBPCs_num%d', ii)), 'png');
    close(fig);

end

%% Save outlines of FDR significance to overlay on consensus spectrograms

xlims = [-0.1, 1];
ttLim = tt(tt>=xlims(1) & tt<xlims(2));

% Apply 5% FDR to visible time range
h_fdr = zeros(size(pSpecs));

%applying FDR across all consensus BPCs at once
pSpecsRelevant = pSpecs(:, tt>=xlims(1) & tt <=xlims(2), :);
h = fdr_bh(pSpecsRelevant(:), 0.05, 'dep', 'yes'); % dependent BH procedure, 1995, FDR = 0.05, apply across all BPCs + exc at same time
h_fdr(:, tt>=xlims(1) & tt <=xlims(2), :) = reshape(h, size(pSpecsRelevant));

for ii = 1:5
    h_fdrCurr = h_fdr(:, tt>=xlims(1) & tt<xlims(2), ii);
    figure; uimagesc(ttLim, f, h_fdrCurr);
    I = getimage(gcf); close(gcf); % get image after interpolation
    I = imresize(I, [size(I, 1)/2, size(I, 2)/4], 'nearest'); % upsample to get tighter boundaries

    mask = boundarymask(I); % get boundaries in even interpolated space
    fig = figure;

    imagesc(ttLim, f, mask);
    colormap([1, 1, 1; 0, 1, 0]);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xlim(xlims);
    ylabel('t-statistic');
    axis xy;
    colorbar;

    if ii == 5
        xlabel('time (s)');
        title('Excluded CCEPs');
    else
        title(sprintf('consensus BPC %d', ii));
    end

    saveas(fig, fullfile('output', 'spectralBB', sprintf('specFDRMask_consBPCs_num%d', ii)), 'svg');
    close(fig);
end

%% Calculate and save time-varying broadband traces for each consensus BPC + excluded group

xlims = [-0.1, 1];
cmapBPC = getCmapVTC('bpc');
cmapBPC(5, :) = [0, 0, 0]; % extra black row for excluded bpcs

BBsCurr = cell(5, 1);
consensusBB = zeros(length(tt), 5); % mean broadband for each consensus BPC
pBBs = zeros(length(tt), 5); % uncorrected t-test pvalue for each consensus broadband

for ii = 1:4
    BBsCurr{ii} = cat(2, masterData{[masterData{:, 3}] == ii & [masterData{:, 4}] >= 1, 7}); % SNR > 1 for each consensus BPC
    consensusBB(:, ii) = log10(geomean(BBsCurr{ii}, 2)); % mean
    [~, p] = ttest(log10(BBsCurr{ii})', 0, 'Alpha', 0.05);
    pBBs(:, ii) = p;
end

% Excluded stim pairs = "BPC #5"
BBsCurr{5} = cat(2, masterData{isnan([masterData{:, 3}]), 7}); % all stim pairs with no consensus BPC label
consensusBB(:, 5) = log10(geomean(BBsCurr{5}, 2)); % mean
[~, p] = ttest(log10(BBsCurr{5})', 0, 'Alpha', 0.05);
pBBs(:, 5) = p;

% Correcting p-values across all consensus + excluded at the same time
h_fdr = zeros(size(pBBs)); % corrected p-values
pBBsAll = pBBs(tt >= xlims(1) & tt <= xlims(2), :); % take chunk across all "5" consensus BPCs in the time interval
h = fdr_bh(pBBsAll(:), 0.05, 'dep', 'yes'); % dependent BH FDR correction
h = reshape(h, size(pBBsAll));
for ii = 1:5
    h_fdr(tt >= xlims(1) & tt <= xlims(2), ii) = h(:, ii);
end

for ii = 1:5
    fprintf('Consensus BPC %d\n', ii);

    fig = figure; hold on
    plotCurvConf(tt, log10(BBsCurr{ii}'), cmapBPC(ii, :), 0.3);
    plot(tt, consensusBB(:, ii), 'LineWidth', 1.5, 'Color', cmapBPC(ii, :));
    h_fdrCurr = h_fdr(:, ii);
    plot(tt(h_fdrCurr' & consensusBB(:, ii)' > 0), 3*h_fdrCurr(h_fdrCurr & consensusBB(:, ii) > 0)', '.', 'MarkerSize', 10, 'Color', 'r');
    plot(tt(h_fdrCurr' & consensusBB(:, ii)' < 0), 3*h_fdrCurr(h_fdrCurr & consensusBB(:, ii) < 0)', '.', 'MarkerSize', 10, 'Color', 'b');
    hold off

    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xlim(xlims); ylim([-1, 4]);
    ylabel('log_{10} broadband power');
    title(sprintf('Consensus BPC %d', ii));

    saveas(fig, fullfile('output', 'spectralBB', sprintf('BB_consBPCs_num%d', ii)), 'svg');
    saveas(fig, fullfile('output', 'spectralBB', sprintf('BB_consBPCs_num%d', ii)), 'png');
    
    hold on % save zoomed in version from 0 to 0.5s
    plot(tt(h_fdrCurr' & consensusBB(:, ii)' > 0), 0.4*h_fdrCurr(h_fdrCurr & consensusBB(:, ii) > 0)', '.', 'MarkerSize', 10, 'Color', 'r');
    plot(tt(h_fdrCurr' & consensusBB(:, ii)' < 0), 0.4*h_fdrCurr(h_fdrCurr & consensusBB(:, ii) < 0)', '.', 'MarkerSize', 10, 'Color', 'b');
    hold off
    xlim([0, 0.5]); ylim([-0.5, 0.5]);

    saveas(fig, fullfile('output', 'spectralBB', sprintf('BBZoomed_consBPCs_num%d', ii)), 'svg');
    saveas(fig, fullfile('output', 'spectralBB', sprintf('BBZoomed_consBPCs_num%d', ii)), 'png');
    close(fig);
end
    

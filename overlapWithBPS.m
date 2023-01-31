%% This script applies BPC method to spectrograms to generate basis profile spectrograms
% Similarity to BPC clusters is quantified by Adjusted Rand Index
% 2022/12/13
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Saves BPSs and compares to consensus BPCs
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
%% Subject names and electrodes

clear;
set(0, 'DefaultFigureRenderer', 'painters');

subNames = {'1', '2', '3', '4', '5'};
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6', 'RC6', 'RC7'}, {'RA9', 'RA10'}};
hemis = 'lrlrr';

%% Load and concatenate voltage BPC outputs and spectrograms

BPCtable = cell(0, 6); % sub, ch, subject BPC, SNR, stimsite, Destrieux label
specData = cell(0, 4);
srate = 2048;

for ii = 1:length(subNames)
    sub = subNames{ii};
    
    for jj = 1:length(subChs{ii})
        ch = subChs{ii}{jj};
        fprintf('sub-%s, ch-%s\n', sub, ch);
        
        % Load spectrogram data from subject
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
        
        % Load voltage BPC SNR data
        snrTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCSNR.tsv', sub, ch)), 'electrical_stimulation_site', 1);  
        stimSites = snrTable.electrical_stimulation_site;
        
        for kk = 1:length(stimSites)
            BPCtable = [BPCtable; {sub}, {ch}, stimSites(kk), snrTable.labelDestrieux(kk), {snrTable.BPC(kk)}, {snrTable.SNR(kk)}];
        end
        
    end
end

BPCtable = cell2table(BPCtable, 'VariableNames', {'sub', 'ch', 'stimSite', 'label', 'voltBPC', 'voltSNR'});

clear S f sigdataBB tt pxx pxxWinds fPSD

%% Calculate basis profile spectrograms

mkdir(fullfile('output', 'BPS'));

rng('default'); % for reproducibility

ttBPC = ttDown(ttDown >= 0.1 & ttDown <= 0.5);
BPCtable.specBPC = nan(height(BPCtable), 1);
BPCtable.specSNR = nan(height(BPCtable), 1);

cmapSpectra = getCmapSpec();

for ii = 1:length(specData)
    events_ch = specData{ii, 3};
    sites = groupby(events_ch.electrical_stimulation_site);
    sites(cellfun(@length, sites(:, 2)) < 8, :) = []; % exclude sites with < 8 trials
    pairTypes = struct('pair', sites(:, 1), 'indices', sites(:, 2));
        
    [B, exc, ~] = bpc_identify(specData{ii, 4}, pairTypes, 1, 50); % 50 reruns, log10 psd
    
    % pull out spectrograms and reshape
    Bspecs = zeros(length(fspec), length(ttBPC), length(B));
    for bb = 1:length(B)
        Bspecs(:, :, bb) = reshape(B(bb).curve, length(fspec), length(ttBPC));
    end
    
    % sort BPCs such that most negative BPC is first, for plotting consistency across electrodes & subjects
    [~, bpcOrd] = sort(squeeze(sum(Bspecs(fspec >= 60 & fspec <= 120, ttBPC >= 0.15 & ttBPC <= 0.35, :), [1, 2])));
    Bspecs = Bspecs(:, :, bpcOrd);
    B = B(bpcOrd);
    
    % add plotweights, reshape each spectrogram and plot
    for bb = 1:length(B)
        
        figure; uimagesc(ttBPC, fspec, Bspecs(:, :, bb), [-0.03, 0.03]);
        colormap(cmapSpectra);
        xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        axis xy;
        colorbar;
        saveas(gcf, fullfile('output', 'BPS', sprintf('sub-%s_ch-%s_specBPC%d', specData{ii, 1}, specData{ii, 2}, bb)), 'png');
        saveas(gcf, fullfile('output', 'BPS', sprintf('sub-%s_ch-%s_specBPC%d', specData{ii, 1}, specData{ii, 2}, bb)), 'svg');
        close(gcf)
       
        B(bb).plotweights = cellfun(@(a, e) mean(a./e.^0.5), B(bb).alphas, B(bb).ep2); % mean alpha/sqrt(ep2) for each group
        for pp = 1:length(B(bb).pairs)
            siteCurr = sites{B(bb).pairs(pp), 1}; % stim site in current BPC group
            
            rowIdx = find(strcmp(BPCtable.sub, specData(ii, 1)) & strcmp(BPCtable.ch, specData(ii, 2)) & strcmp(BPCtable.stimSite, siteCurr));
            assert(length(rowIdx) == 1, 'Found %d rows matching current BPC stim site', length(rowIdx));
            
            BPCtable.specBPC(rowIdx) = bb;
            BPCtable.specSNR(rowIdx) = B(bb).plotweights(pp);
        end
    end
    
end

BPCtable.specBPC(BPCtable.specSNR <= 0) = nan; % exclude small number of stim sites whose SNR <= 0 (no fit)
BPCtable.specSNR(BPCtable.specSNR <= 0) = nan;

% stricter version that only keeps top 50% of global BPS and BPC assignments, based on SNR cutoff
BPCtableStrict = BPCtable;
BPCtableStrict.voltSNR(BPCtable.voltSNR < prctile(BPCtable.voltSNR, 50)) = nan;
BPCtableStrict.voltBPC(BPCtable.voltSNR < prctile(BPCtable.voltSNR, 50)) = nan;
BPCtableStrict.specSNR(BPCtable.specSNR < prctile(BPCtable.specSNR, 50)) = nan;
BPCtableStrict.specBPC(BPCtable.specSNR < prctile(BPCtable.specSNR, 50)) = nan;

%% Determine similarity between voltage BPCs vs BPSs, based on Adjusted Rand Index (ARI)

ARIs = cell(0, 4);
for ii = 1:length(subNames)
    sub = subNames{ii};
    
    for jj = 1:length(subChs{ii})
        ch = subChs{ii}(jj);
        vBPC = BPCtable.voltBPC(strcmp(BPCtable.sub, sub) & strcmp(BPCtable.ch, ch));
        sBPC = BPCtable.specBPC(strcmp(BPCtable.sub, sub) & strcmp(BPCtable.ch, ch));
        vBPC(isnan(vBPC)) = 0; % treat nans as one group for each
        sBPC(isnan(sBPC)) = 0;
        
        ariAll = rand_index(vBPC, sBPC, 'adjusted');
        
        % more selective ARI only including top 80%tile BPCs
        vBPC = BPCtableStrict.voltBPC(strcmp(BPCtableStrict.sub, sub) & strcmp(BPCtableStrict.ch, ch));
        sBPC = BPCtableStrict.specBPC(strcmp(BPCtableStrict.sub, sub) & strcmp(BPCtableStrict.ch, ch));
        vBPC(isnan(vBPC)) = 0; % assign nans as 0-group
        sBPC(isnan(sBPC)) = 0;
        
        ariSig = rand_index(vBPC, sBPC, 'adjusted');
        
        ARIs = [ARIs; {sub}, {ch}, {ariAll}, {ariSig}]; % calculate adjusted rand index        
        
    end
end
ARIs = cell2table(ARIs, 'VariableNames', {'sub', 'ch', 'ARI', 'ARISig'});

% Plot bar graphs of the ARI values
figure('Position', [200, 200, 800, 400]);
bar([ARIs.ARI, ARIs.ARISig], 0.8); ylim([-0.1, 0.8]);
title('ARIs for all subject-electrode pairs');
saveas(gcf, fullfile('output', 'BPS', 'ARIs'), 'svg');
saveas(gcf, fullfile('output', 'BPS', 'ARIs'), 'png');

%% Plot BPS weights on gifti brains

cm = getCmapVTC('bps'); % get colors for BPSs
cm2 = brighten(cm, -0.3);

for ii = 1:length(subNames)
    sub = subNames{ii};
    
    gii = gifti(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', sub), sprintf('pial.%s.surf.gii', upper(hemis(ii)))));
    electrodes = readtableRmHyphens(fullfile('data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg', sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub)));
    xyzs = [electrodes.x, electrodes.y, electrodes.z];
    
    for jj = 1:length(subChs{ii})
        ch = subChs{ii}{jj};
        
        BPCtabCurr = BPCtable(strcmp(BPCtable.sub, sub) & strcmp(BPCtable.ch, ch), :);
        sites = BPCtabCurr.stimSite;
        xyzsPair = ieeg_getPairXyzs(split(sites(:, 1), '-', 2), electrodes);

        f = figure('Position', [1000, 100, 1000, 800]); % normal gifti

        ieeg_RenderGifti(gii); alpha 0.2; hold on
        switch hemis(ii)
            case 'r'
                ieeg_viewLight(90, -40);
            case 'l'
                ieeg_viewLight(-90, -40);
        end

        plot3(xyzs(:,1), xyzs(:,2), xyzs(:,3), 'o', 'Color', 'k', 'MarkerSize', 6, 'MarkerFaceColor', 'w');
        tgt = find(strcmp(electrodes.name, ch)); % circle target electrode
        plot3(xyzs(tgt,1), xyzs(tgt,2), xyzs(tgt,3), 'o', 'Color', 'r', 'MarkerSize', 12, 'LineWidth', 3);
        
        for b = 1:max(BPCtabCurr.specBPC)
            ix_bool = BPCtabCurr.specBPC == b;
            color_add_custom(xyzsPair(ix_bool, :), BPCtabCurr.specSNR(ix_bool), cm2(b, :), 0.6, 1, [8 , 20], 's');
        end
        color_add_custom(xyzsPair(isnan(BPCtabCurr.specBPC), :), 0.5, 0.1*[1 1 1], 0.8, 1, [6, 10], 's'); %non-BPCs
        hold off

        saveas(f, fullfile('output', 'BPS', sprintf('sub-%s_ch-%s_giftispecBPCs', sub, ch)), 'png');
        set(f, 'Position', [1000, 100, 600, 500])
        saveas(f, fullfile('output', 'BPS', sprintf('sub-%s_ch-%s_giftispecBPCs_small', sub, ch)), 'png');
        close(f)
        
    end
end

%% File for calculating consensus BPCs and plotting to MNI brain
% 2021/09/28
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Calculates consensus BPCs from subject-level outputs
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
%% Configure subjects and output dir

clear;
set(0, 'DefaultFigureRenderer', 'painters');

% Subjects 1-4, collateral sulcus measurement electrodes
subNames = {'1', '2', '3', '4'};
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6'}};
hemis = 'lrlr';

mkdir(fullfile('output', 'clusterBPCs'));

%% Load all BPCs into a matrix, X.

X = [];
labels = {};
for ii = 1:length(subNames)
    for jj = 1:length(subChs{ii})
        curvesObj = load(fullfile('output', sprintf('sub-%s', subNames{ii}), sprintf('sub-%s_ch-%s_BPCs.mat', subNames{ii}, subChs{ii}{jj})));
        X = [X; curvesObj.curves'];
        for kk = 1:size(curvesObj.curves, 2)
            labels = [labels; sprintf('%s-%s-%d', subNames{ii}, subChs{ii}{jj}, kk)];
        end
    end
end

tt = curvesObj.ttBPC;

%% save the matrix of all BPCs across subjects

Xunweight = weightExponential(X', tt, 0.1, true);
plotTrials(tt, Xunweight, 0.1, [], [], 'LineWidth', 1);
xlabel('Time from Stimulation (s)'); ylabel('Voltage (Unit-Normalized)');

saveas(gcf, fullfile('output', 'clusterBPCs', 'allBPCs'), 'svg');

%% SVD (non-centered PCA) on X

[loadings, score, ~, ~, explained] = pca(X, 'Centered', false); % no centering because baseline already corrected

% Plot variance explained
figure; plot(cumsum(explained), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k'); xlabel('num PCs'); ylabel('% variance explained'); 
xlim([0, 20]);
ylim([0, 100]);
saveas(gcf, fullfile('output', 'clusterBPCs', 'varExp'), 'png');

% plot selected PC loadings
figure;
toPlot = weightExponential(loadings(:, 1:3), tt, 0.1, true)';
ys = plotTrials(tt, toPlot', 0.2, [], repmat([0 0 0], [3, 1]), 'LineWidth', 1.5);
xlabel('time (s)'); ylabel('PC loadings');
set(gca, 'YTick', flip(ys)); set(gca, 'YTicklabels', {'PC3', 'PC2', 'PC1'});
saveas(gcf, fullfile('output', 'clusterBPCs', 'loadings'), 'png');

% plot BPCs in PC space.
figure('Position', [1000, 600, 800, 600]);
plot3(score(:, 1), score(:, 2), score(:, 3), 'k.', 'MarkerSize', 18);
view(50, 30);
text(score(:, 1), score(:, 2), score(:, 3), labels); % labels are "sub-electrode-BPC#"
axis equal
xlabel('SV1'); ylabel('SV2'); zlabel('SV3');
saveas(gcf, fullfile('output', 'clusterBPCs', 'PCspace'), 'png');

%% K-means clustering in SV space

cm = getCmapVTC('bpc');
k = 4; % number of clusters
nPCs = 3; % number of PCs to cluster on
[idx, C] = kmeans(score(:, 1:nPCs), k, 'Replicates', 100); % get kmeans clusters

% sort clusters by centroid position along PC1
[~, order] = sort(C(:, 1));
idx = changem(idx, 1:k, order);
C = C(order, :);

% Plot centroid linear combs of loadings as the consensus BPCs
Vcent = loadings(:, 1:nPCs)*C';
Vcent = weightExponential(Vcent, tt, 0.1, true); % antiweight for plotting
figure;
ys = plotTrials(tt, Vcent, 0.2, [], cm, 'LineWidth', 1.5);
xlabel('Time from Stimulation (s)'); ylabel('Voltage (Unit-Normalized)');
saveas(gcf, fullfile('output', 'clusterBPCs', 'consensusBPCs'), 'png');
saveas(gcf, fullfile('output', 'clusterBPCs', 'consensusBPCs'), 'svg');

% Plot BPCs in PC space, labeled by consensus BPC color
figure('Position', [1000, 600, 800, 600]); hold on
for i=1:size(C, 1)
    plot3(score(idx==i, 1), score(idx==i, 2), score(idx==i, 3), 'o', 'MarkerFaceColor', cm(i, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
end
xticks([]); yticks([]); zticks([]);
xlabel('SV1'); ylabel('SV2'); zlabel('SV3');
view(50, 30);
axis equal
hold off
saveas(gcf, fullfile('output', 'clusterBPCs', 'PCspaceLabelled'), 'png');
saveas(gcf, fullfile('output', 'clusterBPCs', 'PCspaceLabelled'), 'svg');

% Find peaks of each consensus BPC, print out
[~, locs] = findpeaks(-Vcent(:, 1), 'MinPeakProminence', 0.01);
fprintf('Consensus BPC1: neg peak at %.0f ms\n', 1000*tt(locs(1)));
[~, samp1] = min(Vcent(:, 2)); [~, locs] = findpeaks(Vcent(:, 2), 'MinPeakProminence', 0.01);
fprintf('Consensus BPC2: neg peak at %.0f ms, pos peak at %.0f ms, pos peak at %.0f ms\n', 1000*tt(samp1), 1000*tt(locs(end-1)), 1000*tt(locs(end)));
[~, samp1] = max(Vcent(:, 3)); [~, samp2] = min(Vcent(:, 3)); [~, samp3] = findpeaks(Vcent(:, 3), 'MinPeakProminence', 0.01); 
fprintf('Consensus BPC3: pos peak at %.0f ms, neg peak at %.0f ms, pos peak at %.0f ms\n', 1000*tt(samp1), 1000*tt(samp2), 1000*tt(samp3));
[~, samp] = max(Vcent(:, 4));
fprintf('Consensus BPC4: pos peak at %.0f ms\n', 1000*tt(samp));

%% Plot all stim sites to MNI brain, colored by consensus BPC
% - Compute anatomical distribution bar graphs
% - Save consensus BPC label as column to original BPCSnr file

% load MNI 152 brain and put into Gifti format
load(fullfile('data', 'derivatives', 'MNI_cortex_left.mat'));
cortexGii = struct();
cortexGii.vertices = cortex.vert; cortexGii.faces = cortex.tri;
clear cortex;

figure('Position', [1000, 100, 700, 500]);
ieeg_RenderGifti(cortexGii); ieeg_viewLight(-90, -45); alpha 0.3; hold on

countBPCs = 1; % counter for BPCs across subjects
allConsensusLabels = []; % used to cross-tabulate stim region by consensus BPC
allStimRegions = [];
masterMNItab = {}; % sub, elec, stimsite, x, y, z, consBPC, SNR
for ii = 1:length(subNames)
    elecsMNI = readtableRmHyphens(fullfile('output', sprintf('sub-%s', subNames{ii}), 'MNI152', sprintf('sub-%s_ses-ieeg01_space-MNI152NLin6Sym_electrodes.tsv', subNames{ii}))); % created with saveMNIelectrodes.m
    if strcmp(hemis(ii), 'r'), elecsMNI.x = -elecsMNI.x; end % flip for right hemis so as to plot on left hemisphere
    
    for jj = 1:length(subChs{ii})
        snrTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', subNames{ii}), sprintf('sub-%s_ch-%s_BPCSNR.tsv', subNames{ii}, subChs{ii}{jj})), 'electrical_stimulation_site', 1);
        recElec = elecsMNI(strcmp(elecsMNI.name, subChs{ii}{jj}), :); % recording electrode
        pairXyzs = ieeg_getPairXyzs(split(snrTable.electrical_stimulation_site, '-', 2), elecsMNI);
        
        % replace subject-level BPC labels with consensus BPC labels
        bpcLabels = snrTable.BPC(~isnan(snrTable.BPC));
        SNR = snrTable.SNR(~isnan(snrTable.SNR));
        destrieux = snrTable.labelDestrieux(~isnan(snrTable.SNR));
        consensusBpcLabels = changem(bpcLabels, idx(countBPCs:countBPCs+length(unique(bpcLabels))-1), unique(bpcLabels));
        
        % save consensus BPC labels to BPCSNR file
        snrTable(:, 5:end) = [];
        snrTable.consBPC = nan(height(snrTable), 1);
        snrTable.consBPC(~isnan(snrTable.SNR)) = consensusBpcLabels;
        writetable(snrTable, fullfile('output', sprintf('sub-%s', subNames{ii}), sprintf('sub-%s_ch-%s_BPCSNR.tsv', subNames{ii}, subChs{ii}{jj})), 'Filetype','text', 'Delimiter','\t');
                            
        allConsensusLabels = [allConsensusLabels; consensusBpcLabels(SNR >= 1)];
        allStimRegions = [allStimRegions; getStimRegion(destrieux(SNR >= 1))];
        for kk = 1:height(snrTable)
            masterMNItab = [masterMNItab;
                {subNames{ii}, subChs{ii}{jj}, snrTable.electrical_stimulation_site(kk), pairXyzs(kk, 1), pairXyzs(kk, 2), pairXyzs(kk, 3), snrTable.consBPC(kk), snrTable.SNR(kk)}]; % all stim sites
        end
        
        countBPCs = countBPCs + length(unique(bpcLabels)); % update BPC counter
        
        plot3(recElec.x, recElec.y, recElec.z, 'o', 'Color', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        uniqueConsensusLabels = unique(consensusBpcLabels);
        for kk = 1:length(uniqueConsensusLabels)
            color_add_custom(pairXyzs(consensusBpcLabels == uniqueConsensusLabels(kk) & SNR >= 1, :), ...
                             SNR(consensusBpcLabels == uniqueConsensusLabels(kk) & SNR >= 1), cm(uniqueConsensusLabels(kk), :), 0.8, 3, [6, 20], 's');
        end
        color_add_custom(pairXyzs(isnan(snrTable.BPC) | snrTable.SNR < 1, :), 0.5, 0.1*[1, 1, 1], 0.8, 1, [1, 10], 's');
        
    end
end
saveas(gcf, fullfile('output', 'clusterBPCs', 'giftiConsBPCs'), 'png');

% Create distribution bar graphs and display the contingency table consensus BPC category x anatomic bin
[tbl, chi2, p_cont] = crosstab(allConsensusLabels, allStimRegions); % contingency table, chi2 independence test
disp(tbl); fprintf('chi-squared P = %.2e\n', p_cont);

cmAnat = getCmapVTC('anat');
f = figure('Position', [800, 600, 600, 600]);
b = bar(categorical(unique(allConsensusLabels)), tbl, 'stacked', 'BarWidth', 0.7);
for ii = 1:length(b)
     b(ii).FaceColor = cmAnat(ii, :);
end
ylim([0, max(sum(tbl, 2)) + 1]);
xlabel('Consensus BPC'); ylabel('number of stim sites');
title('Consensus BPC distribution all subjects');
saveas(gcf, fullfile('output', 'clusterBPCs', 'barGraphsConsBPC'), 'svg');
saveas(gcf, fullfile('output', 'clusterBPCs', 'barGraphsConsBPC'), 'png');

% print number of stim sites without BPC
fprintf('%d stim sites with no BPC category\n', sum(isnan([masterMNItab{:, 8}])));

% print number of excluded stim sites based on SNR < 1
fprintf('%d stim sites out of %d excluded by SNR = 1 threshold\n', sum([masterMNItab{:, 8}] < 1), sum(~isnan([masterMNItab{:, 8}])));

% single tsv that consolidates all stim site coordinates in MNI 152 space and consensus BPC categories
writetable(cell2table(masterMNItab, 'VariableNames', {'sub', 'elec', 'stimsite', 'x', 'y', 'z', 'consBPC', 'SNR'}), ...
    fullfile('output', 'clusterBPCs', 'masterMNItab.tsv'), 'Filetype','text', 'Delimiter','\t');



%% For figure 3-1: Leave one subject out cross-validation analysis.

mkdir(fullfile('output', 'clusterBPCs', 'xval'));

% Reload BPCs by subject
Xcross = cell(length(subNames), 1);
idxTrue = cell(length(subNames), 1); % "true" consensus BPC labels, separated by subjects, from original non cross validated model
idxTest = cell(length(subNames), 1); % labels from cross validation

cc = 1;
for ii = 1:length(subNames)
    for jj = 1:length(subChs{ii})
        curvesObj = load(fullfile('output', sprintf('sub-%s', subNames{ii}), sprintf('sub-%s_ch-%s_BPCs.mat', subNames{ii}, subChs{ii}{jj})));
        Xcross{ii} = [Xcross{ii}; curvesObj.curves'];
    end
    
    idxTrue{ii} = idx(cc:(cc+size(Xcross{ii}, 1)-1)); % set of all true consensus BPC assignments belonging to current subject
    cc = cc + size(Xcross{ii}, 1);
end

% Concatenate 3/4 subjects to make training set on each CV fold
for ii = 1:length(subNames)
    
    fprintf('Test subject = sub-%s\n', subNames{ii});
    
    boolTrain = ~antifind(ii, length(subNames));
    Xtest = Xcross{ii}; % BPCs for test subject
    Xtrain = cat(1, Xcross{boolTrain}); % BPCs for training subjects
    
    % save plot of which BPCs were used for test then train set
    figure; plotTrials(tt, weightExponential(Xtest', tt, 0.1, true), 0.1, [], [], 'LineWidth', 1);
    xlabel('Time from Stimulation (s)'); ylabel('Voltage (Unit-Normalized)');
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('BPCsTest%d', ii)), 'svg');
    figure; plotTrials(tt, weightExponential(Xtrain', tt, 0.1, true), 0.1, [], [], 'LineWidth', 1);
    xlabel('Time from Stimulation (s)'); ylabel('Voltage (Unit-Normalized)');
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('BPCsTrain%d', ii)), 'svg');
    
    % SVD on training set
    [loadingsTrain, scoreTrain] = pca(Xtrain, 'Centered', false);

    % Kmeans
    k = 4; % number of clusters
    nPCs = 3; % number of PCs to cluster on
    [idxTrain, CTrain] = kmeans(scoreTrain(:, 1:nPCs), k, 'Replicates', 100); % get kmeans clusters (group assignments and centroids)

    % sort clusters by centroid order (similarity) to full model consensus BPC 1, so that the labels present in consistent order regardless of PCAs formed
    VcentTrain = loadingsTrain(:, 1:nPCs)*CTrain';
    sim = dot(VcentTrain, repmat(loadings(:, 1), 1, k)); % similarity to the first PC from the overall model
    [~, order] = sort(sim); % assign consensus labels in ascending similarity
    idxTrain = changem(idxTrain, 1:k, order);
    CTrain = CTrain(order, :);
    VcentTrain = VcentTrain(:, order);

    % Plot centroid linear combs of loadings == consensus BPC shapes
    VcentTrain = weightExponential(VcentTrain, tt, 0.1, true); % antiweight for plotting
    figure;
    ys = plotTrials(tt, VcentTrain, 0.15, [], cm, 'LineWidth', 1);
    xlabel('Time from Stimulation (s)'); ylabel('Voltage (Unit-Normalized)');
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('consBPCsTrain%d', ii)), 'png');
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('consBPCsTrain%d', ii)), 'svg');
    
    % Plot training BPCs in PC space 
    figure('Position', [1000, 600, 800, 600]); hold on
    for jj = 1:size(CTrain, 1)
        plot3(scoreTrain(idxTrain == jj, 1), scoreTrain(idxTrain == jj, 2), scoreTrain(idxTrain == jj, 3), 'o', 'MarkerFaceColor', cm(jj, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
    end
    xticks([]); yticks([]); zticks([]);
    xlabel('SV1'); ylabel('SV2'); zlabel('SV3');
    view(50, 40);
    axis equal
    hold off
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('PCspaceLabTrain%d', ii)), 'png');
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('PCspaceLabTrain%d', ii)), 'svg');
    
    % Replot with test BPCs labeled
    figure('Position', [1000, 600, 800, 600]); hold on
    for jj = 1:size(CTrain, 1)
        plot3(scoreTrain(idxTrain == jj, 1), scoreTrain(idxTrain == jj, 2), scoreTrain(idxTrain == jj, 3), 'o', 'MarkerFaceColor', cm(jj, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
    end
    xticks([]); yticks([]); zticks([]);
    xlabel('SV1'); ylabel('SV2'); zlabel('SV3');
    view(50, 40);
    
    % get labels for test set
    PCtest = Xtest*loadingsTrain(:, 1:nPCs); % PC space locations for testing set    
    for jj = 1:size(PCtest, 1)
        dists = vecnorm(PCtest(jj, :) - CTrain, 2, 2); % euclidean distance to each training centroid
        [~, idxTest{ii}(jj)] = min(dists);
        plot3(PCtest(jj, 1), PCtest(jj, 2), PCtest(jj, 3), 's', 'MarkerFaceColor', cm(idxTest{ii}(jj), :), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
    end
    idxTest{ii} = idxTest{ii}(:); % ensure column vector
    
    axis equal
    hold off
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('PCspaceLabTrainTest%d', ii)), 'png');
    saveas(gcf, fullfile('output', 'clusterBPCs', 'xval', sprintf('PCspaceLabTrainTest%d', ii)), 'svg');
    
end

% Assess per-subject test set accuracy
acc = zeros(length(subNames), 1);
for ii = 1:length(subNames)
    acc(ii) = sum(idxTest{ii} == idxTrue{ii})/length(idxTrue{ii});
    fprintf('Accuracy on test sub %s = %0.1f%%\n', subNames{ii}, acc(ii)*100);
end
accTotal = sum(cat(1, idxTest{:}) == idx)/length(idx);

fprintf('Total accuracy = %0.1f%%\n', accTotal*100);



%% Functions

% Used to group stim site regions by destrieux label
function stimRegions = getStimRegion(destrieuxLabels)
    stimRegions = 6*ones(length(destrieuxLabels), 1); % "other" as number 6 so it comes last
    stimRegions(contains(lower(destrieuxLabels), 'hippocampus')) = 1;
    stimRegions(contains(lower(destrieuxLabels), 'amygdala')) = 2;
    stimRegions(contains(lower(destrieuxLabels), {'oc-temp', 'collat'})) = 3;
    stimRegions(contains(lower(destrieuxLabels), {'temporal_', 'temp_sup'})) = 4;
    stimRegions(contains(lower(destrieuxLabels), 'ins')) = 5;
end

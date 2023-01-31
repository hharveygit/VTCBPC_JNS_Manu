%% File examining relationship between BPC SNR and distance to measurement electrode
% 2023/01/11
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package
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
%% Configure subjects (1-4, consensus BPC electrodes)

clear;
set(0, 'DefaultFigureRenderer', 'painters');

% consensus BPC electrodes
subNames = {'1', '2', '3', '4'};
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6'}};
hemis = 'lrlr';

subs = struct();
for ii = 1:length(subNames)
    subs(ii).name = subNames{ii};
    subs(ii).chs = subChs{ii};
    subs(ii).hemi = hemis(ii);
end

mkdir(fullfile('output', 'distanceVBPC'));

clear subNames hemis subChs

%% Calculate distance between measurement electrode and stim site for all collateral sulcus electrode stim sites

% Calculate MNI distance from measurement electrode across all subjects and electrodes
masterMniTab = readtable(fullfile('output', 'clusterBPCs', 'masterMNItab.tsv'), ....
                         'FileType', 'text', 'Delimiter', '\t'); % has all the MNI position, SNR, consensus BPC labels across subjects and electrodes
                     
xyzsStimMni = [masterMniTab.x, masterMniTab.y, masterMniTab.z];
distance = nan(height(masterMniTab), 1); % MNI distance from recording electrode

elecsXyzMni = []; % store the MNI coordinate of all measurement electrodes to plot in next section
for ii = 1:length(subs)
    
    sub = subs(ii).name;

    % Load MNI 152 electrode positions, reflect x for right hemisphere electrodes to match the reflected values stored for stim sites in masterMniTab
    electrodesMni = readtableRmHyphens(fullfile('output', sprintf('sub-%s', sub), 'MNI152', sprintf('sub-%s_ses-ieeg01_space-MNI152NLin6Sym_electrodes.tsv', sub)));
    if strcmp(subs(ii).hemi, 'r'), electrodesMni.x = -electrodesMni.x; end
    
    xyzsMni = [electrodesMni.x, electrodesMni.y, electrodesMni.z];
    
    for jj = 1:length(subs(ii).chs)
        
        ch = subs(ii).chs{jj};
        chXyzMni = xyzsMni(strcmp(electrodesMni.name, ch), :); % measurement electrode position
        elecsXyzMni = [elecsXyzMni; chXyzMni];
        
        stimSitesCurr = masterMniTab.sub == str2double(sub) & strcmp(masterMniTab.elec, ch); % bool indices of stim site rows belonging to this channel
        distance(stimSitesCurr) = vecnorm(xyzsStimMni(stimSitesCurr, :) - chXyzMni, 2, 2); % euclidean MNI distance        
    end
end

assert(~any(isnan(distance)), 'Distances not calculated for some subjects/electrodes');
masterMniTab.distance = distance;

%% 1. Plot distance from measurement elec to stim site, grouped by consensus BPC. Q: Is distance different across consensus BPC groups?

% Plot BPC assignment vs distance, test group mean differences with 1-way ANOVA
bpclabel = masterMniTab.consBPC; bpclabel(isnan(bpclabel)) = 5; % consensus BPC as grouping variable, put nan at the end
[p, t, stats] = anova1(masterMniTab.distance, bpclabel, 'off'); % 1-way anova on distance across BPCs
[c, m] = multcompare(stats, 'Display', 'off'); % post-hoc test of means
cTab = array2table(c, "VariableNames", ["group1", "group2", "Lower Limit", "group1-group2", "Upper Limit", "P-value"])
mTab = array2table(m, "VariableNames", ["mean", "Serr"])
writetable(cTab, fullfile('output', 'distanceVBPC', 'consBPCVMniDistancePostHoc.tsv'), 'FileType', 'text', 'Delimiter', '\t');

cm = getCmapVTC('bpc');
cmCurr = [cm; [0.3, 0.3, 0.3]]; % custom colormap to match all consensus BPCs + excluded

figure('Position', [200, 200, 300, 500]);
distributionPlot(masterMniTab.distance, 'groups', bpclabel, 'color', num2cell(cmCurr, 2), 'showMM', 0);
hold on
plotSpread(masterMniTab.distance, 'distributionIdx', bpclabel, 'xValues', 1:5, 'spreadWidth', 1, 'distributionColors', [0.8, 0.8, 0.8]);
xlim([0, 6]); ylim([0, inf]);
title('Consensus BPCs');
text(0.5, max(masterMniTab.distance)+5, sprintf('P=%0.2e', p));
hold off
saveas(gcf, fullfile('output', 'distanceVBPC', 'consBPCVMniDistance'), 'png');
saveas(gcf, fullfile('output', 'distanceVBPC', 'consBPCVMniDistance'), 'svg');

%% 2. Plot distance vs. SNR. Q: Within each consensus BPC, and across all consensus BPCs, is distance correlated with SNR?
% Mixed effects model to model subject as random effect

% scatterplot of SNR vs distance, colored by BPC
x = masterMniTab.distance(~isnan(masterMniTab.consBPC));
y = masterMniTab.SNR(~isnan(masterMniTab.consBPC));
bpclabel = masterMniTab.consBPC(~isnan(masterMniTab.consBPC));
fid = fopen(fullfile('output', 'distanceVBPC', 'consBPC_snrDistCorr.txt'), 'w');
fprintf(fid, 'BPC\tn\tR\tP\n');
[r, p] = corr(x, y); % correlation across all BPCs of distance vs SNR
fprintf(fid, 'All\t%d\t%0.2f\t%0.2e\n', length(x), r, p); % overall correlation
figure('Position', [400, 200, 600, 400]); hold on
for bb = 1:4
    plot(x(bpclabel == bb), y(bpclabel == bb), '.', 'Color', cm(bb, :), 'MarkerSize', 14);
    [r, p] = corr(x(bpclabel == bb), y(bpclabel == bb)); % correlation across all BPCs of distance vs SNR
    fprintf(fid, '%d\t%d\t%0.2f\t%0.2e\n', bb, sum(bpclabel == bb), r, p); % correlation for current consensus BPC
end
fclose(fid);
saveas(gcf, fullfile('output', 'distanceVBPC', 'consBPCSnrDistance'), 'png');
saveas(gcf, fullfile('output', 'distanceVBPC', 'consBPCSnrDistance'), 'svg');
close all

% Fit a mixed effects model for each cons BPC, with random effect of subject on slope and intercept (and correlation) considered
lmeTab = masterMniTab; lmeTab(isnan(masterMniTab.consBPC), :) = []; % copy table and remove nans, for ease
for bb = 1:4
    lmeTabCurr = lmeTab(lmeTab.consBPC == bb, :);
    lme = fitlme(lmeTabCurr, 'SNR~distance+(distance|sub)');
    fprintf('Consensus BPC %d:\n', bb);
    display(lme.Coefficients);
end

%% Plot all stim sites on MNI brain, not thresholding by minimum SNR threshold

load(fullfile('data', 'derivatives', 'MNI_cortex_left.mat'));
cortexGii = struct();
cortexGii.vertices = cortex.vert; cortexGii.faces = cortex.tri;
clear cortex;

figure('Position', [1000, 100, 700, 500]);
ieeg_RenderGifti(cortexGii); ieeg_viewLight(-90, -45); alpha 0.3; hold on

for bb = 1:4 % add all consensus BPC stim site labels without thresholding
    color_add_custom(xyzsStimMni(masterMniTab.consBPC == bb, :), masterMniTab.SNR(masterMniTab.consBPC == bb), cm(bb, :), 0.8, 3, [6, 20], 's');
end
color_add_custom(xyzsStimMni(isnan(masterMniTab.consBPC), :), 0.5, 0.1*[1, 1, 1], 0.8, 1, [1, 10], 's'); % excluded bpcs

for ii = 1:size(elecsXyzMni, 1)
    plot3(elecsXyzMni(ii, 1), elecsXyzMni(ii, 2), elecsXyzMni(ii, 3), 'o', 'Color', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

hold off
saveas(gcf, fullfile('output', 'distanceVBPC', 'giftiConsBPCsNoThreshold'), 'png');
    
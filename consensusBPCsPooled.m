%% File for calculating consensus BPCs and plotting to MNI brain
%
% 2023/04
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
%% Configure user parameters, paths

set(0, 'DefaultFigureRenderer', 'painters');

subNames = {'1', '2', '3', '4'}; % does not include subject 7 (testing set)
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6'}};
hemis = 'lrlr';

subs = struct();
for ii = 1:length(subNames)
    subs(ii).name = subNames{ii};
    subs(ii).chs = subChs{ii};
    subs(ii).hemi = hemis(ii);
end

mkdir(fullfile('output', 'consensusBPCsPooled'));

% Load data if this part alreadyrun
%load('output/consensusBPCspooled/CCEPsPooled.mat');

%% Load voltage data from all subjects (if not saved)
% divide by SD at each electrode, then pool

% pooled "events_ch" across all subject-electrode combos
eventsPooled = cell2table(cell(0, 4), 'VariableNames', {'subject', 'electrode', 'electrical_stimulation_site', 'electrical_stimulation_current'});

sigdataPooled = [];

for ii = 1:length(subNames)
    sub = subNames{ii};
    chs = subs(ii).chs;
    
    fprintf('sub-%s\n', sub);
    subDir = fullfile(pwd, 'data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg');
    
    % this section copied from saveBPCs.m
    mefPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_ieeg.mefd', sub));
    channelsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_channels.tsv', sub));
    eventsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_events.tsv', sub));

    mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
    mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'});
    mefObj.loadMefTrials([-1, 1.5]); % change to end at 1.5 s for tighter window to normalize around stim onset
    mefObj.car(true); % by 64-block
    mefObj.pruneChannels(chs);
    mefObj.removeLN('SpectrumEstimation');
    mefObj.subtractBaseline([-0.5, -0.05], 'median');

    events = mefObj.evts;
    tt = mefObj.tt;
    srate = mefObj.srate;
    data = mefObj.data;

    % Remove all events that don't fit 'eln-el(n+1)'
    pairNum = zeros(size(events, 1), 2);
    for kk = 1:length(pairNum)
        pairNum(kk, :) = str2double(regexp(events.electrical_stimulation_site{kk}, '\d*', 'match'));
    end
    data(:, :, diff(pairNum, 1, 2) ~= 1) = [];
    events(diff(pairNum, 1, 2) ~= 1, :) = [];

    % Remove all events that correspond to seizure onset zones
    electrodes = readtableRmHyphens(fullfile(subDir, sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub)));
    elecsSoz = electrodes.name(contains(electrodes.seizure_zone, 'SOZ'));
    eventsSoz = any(ismember(split(events.electrical_stimulation_site, '-'), elecsSoz), 2);
    fprintf('Removed %d events at %d sites in SOZ\n', sum(eventsSoz), length(unique(events.electrical_stimulation_site(eventsSoz))));
    data(:, :, eventsSoz) = [];
    events(eventsSoz, :) = [];

    assert(size(data, 3) == height(events), 'data - events mismatch');
    
    for jj = 1:length(subChs{ii})
        
        %% Get data from channel, exclude stim trials at hannel
        ch = chs{jj};
        
        fprintf('\tch-%s\n', ch);
        
        sigdata = squeeze(data(strcmp(mefObj.channels.name, ch), :, :))';
        events_ch = events;
    
        % remove stimulated channels from current events
        sigdata(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
        events_ch(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
        sites = groupby(events_ch.electrical_stimulation_site);

        % resolve only 4mA or 6mA trials for stim sites with both
        for kk = 1:size(sites, 1)
            idxes = sites{ii, 2};
            trialsCurr = events_ch(idxes, :);
            if length(unique(trialsCurr.electrical_stimulation_current)) == 1, continue; end

            trials6ma = idxes(strcmp(trialsCurr.electrical_stimulation_current, '6.0 mA')); % trials at 6 mA stim
            if length(trials6ma) >= 8
                sites{kk, 2} = trials6ma;
                continue
            else % assign the 4.0 mA trials
                sites{kk, 2} = idxes(strcmp(trialsCurr.electrical_stimulation_current, '4.0 mA'));
            end
        end
        sites(cellfun(@length, sites(:, 2)) < 8, :) = []; % exclude sites with < 8 trials 
        
        % propagate only the trials to keep to events_ch. Necessary here because we are concatenating the events_ch and sigdata variables
        idxes2Keep = sort(vertcat(sites{:, 2}));
        events_ch = events_ch(idxes2Keep, :);
        sigdata = sigdata(idxes2Keep, :);
        
        % concatenate the events info about subjects and electrodes
        eventsPooled = [eventsPooled; [repmat({sub}, height(events_ch), 1), repmat({ch}, height(events_ch), 1), events_ch.electrical_stimulation_site, events_ch.electrical_stimulation_current]];
        
        % normalize by SD and concatenate. Note that this is normalizing across entire data range ([-1, 1.5]), which should be a better estimate of baseline.
        sigdata = sigdata / std(sigdata, 0, 'all');
        sigdataPooled = [sigdataPooled; sigdata];
        
    end
end

% Save the pooled data for direct loading
save('output/consensusBPCspooled/CCEPsPooled.mat', 'sigdataPooled', 'tt', 'eventsPooled');

%% Calculate BPCs

rng('default');
%
% Only join sub-stimsite, which means all electrodes in the same subject are considered as one "super-electrode" that should in theory
%       record the same data when an electrode pair is stimulated. This approach is more consistent with the intended use of the BPC algorithm,
%       but doesn't match behavior of original consensus BPCs.
sites = groupby(join([eventsPooled.subject, eventsPooled.electrical_stimulation_site], '-', 2));

disp('Getting BPCs');
tau = 0.1;

pairTypes = struct('pair', sites(:, 1), 'indices', sites(:, 2));
V = sigdataPooled(:, tt >= 0.011 & tt < 0.5)';
ttBPC = tt(tt >= 0.011 & tt < 0.5);
V = weightExponential(V, tt(tt >= 0.011 & tt < 0.5), tau); % weight by decaying exponential of time constant 100ms

[B, exc, H] = bpc_identify(V, pairTypes, 1.2, 50); % 50 reruns

curves = [B.curve];

% Sort BPCs by peak voltage before 50 ms (negative first)
[~, bOrd] = sort(min(curves));
B = B(bOrd);
curves = curves(:, bOrd);

BPCs_expanded = nan(size(sites)); % BPC label and plotweight for each site in same order as sites
for b = 1:length(B)
    B(b).plotweights = cellfun(@(a, e) mean(a./e.^0.5), B(b).alphas, B(b).ep2); % mean alpha/sqrt(ep2) for each group
    BPCs_expanded(B(b).pairs, 1) = b;
    BPCs_expanded(B(b).pairs, 2) = B(b).plotweights;
end

% plot BPCs and save svg
f = figure('Position',[100 100 800 600]); hold on
xlabel('time (s)');
ylabel('V (unit-normalized)');
set(gca, 'YTicklabel', []);
cm = getCmapVTC('bpc');
plotTrials(ttBPC, weightExponential(curves, ttBPC, tau, true), 0.2, [], cm, 'LineWidth', 1.5);
saveas(gcf, fullfile('output', 'consensusBPCspooled', 'consensusBPCs'), 'png');
saveas(gcf, fullfile('output', 'consensusBPCspooled', 'consensusBPCs'), 'svg');

close(f)


%% Plot all stim sites to MNI brain, colored by consensus BPC
% - Compute anatomical distribution bar graphs

load(fullfile('data', 'derivatives', 'MNI_cortex_left.mat'));
cortexGii = struct();
cortexGii.vertices = cortex.vert; cortexGii.faces = cortex.tri;
clear cortex;

% calculate the MNI coordinates for all stim pairs and measurement electrodes
recElecMni = [];
pairXyzsMni = nan(size(sites, 1), 3); % MNI xyz coordinate for each stim site (in all subjects)
for ii = 1:length(subNames)
    sub = subNames{ii};
        
    elecsMNI = readtableRmHyphens(fullfile('output', sprintf('sub-%s', sub), 'MNI152', sprintf('sub-%s_ses-ieeg01_space-MNI152NLin6Sym_electrodes.tsv', sub))); % created with saveMNIelectrodes.m
    if strcmp(hemis(ii), 'r'), elecsMNI.x = -elecsMNI.x; end % flip for right hemis so as to plot on left hemisphere
    
    sitesThisSub = contains(sites(:, 1), sub); % rows of sites from current subject
    pairXyzsMni(sitesThisSub, :) = ieeg_getPairXyzs(extractAfter(sites(sitesThisSub, 1), '-'), elecsMNI); % get stim pair MNI coordinates for all stim pairs in current subject
    
    for jj = 1: length(subChs{ii})
        recElecMni = [recElecMni; elecsMNI(strcmp(elecsMNI.name, subChs{ii}{jj}), :)]; % concatenate recording electrode
    end
    
end
assert(all(~isnan(pairXyzsMni), 'all'), 'Error: missing values not filled in pairXyzsMni');

figure('Position', [1000, 100, 700, 500]);
ieeg_RenderGifti(cortexGii); ieeg_viewLight(-90, -45); alpha 0.3; hold on

plot3(recElecMni.x, recElecMni.y, recElecMni.z, 'o', 'Color', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

for kk = 1:length(B)
    color_add_custom(pairXyzsMni(BPCs_expanded(:, 1) == kk, :), ...
                     BPCs_expanded(BPCs_expanded(:, 1) == kk, 2), cm(kk, :), 0.8, 3, [6, 20], 's');
end
color_add_custom(pairXyzsMni(isnan(BPCs_expanded(:, 1)), :), 0.5, 0.1*[1, 1, 1], 0.8, 1, [1, 10], 's');
saveas(gcf, fullfile('output', 'consensusBPCspooled', 'giftiConsBPCs'), 'png');


%% Quantify anatomical distribution of consensus BPCs

% find destrieux atlas for each stim sites
anatLabels = cell(size(sites, 1), 1);
for ii = 1:size(sites, 1)
    sub = extractBefore(sites{ii, 1}, '-');
    stimsite = extractAfter(sites{ii, 1}, '-');
    
    chs = subChs{strcmp(subNames, sub)};
    
    for jj = 1:length(chs) % try searching through all recording channels to find the stim site
        %fprintf('Reading SNR table from sub-%s, ch-%s\n', sub, chs{jj})
        snrTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCSNR.tsv', sub, chs{jj})), 'electrical_stimulation_site', 1);
        lab = snrTable.labelDestrieux(strcmp(snrTable.electrical_stimulation_site, stimsite));
        
        % keep searching through recording electrodes until stim site found
        if isempty(lab), continue;
        else
            anatLabels(ii) = lab;
            break
        end
        
    end
end

% determine stim regions from the anatomic labels
stimRegions = getStimRegion(anatLabels);

% contingency table and chi2 test of independence
[tbl, chi2, p_cont] = crosstab(BPCs_expanded(BPCs_expanded(:, 2) >= 1, 1), stimRegions(BPCs_expanded(:, 2) >= 1)); % apply only to those with SNR >= 1
disp(tbl); fprintf('chi-squared P = %.2e\n', p_cont);

cmAnat = getCmapVTC('anat');
f = figure('Position', [800, 600, 600, 600]);
b = bar(categorical(1:length(B)), tbl, 'stacked', 'BarWidth', 0.7);
for ii = 1:length(b)
     b(ii).FaceColor = cmAnat(ii, :);
end
ylim([0, max(sum(tbl, 2)) + 1]);
xlabel('Consensus BPC'); ylabel('number of stim sites');
title('Consensus BPC distribution all subjects');
saveas(gcf, fullfile('output', 'consensusBPCsPooled', 'barGraphsConsBPC'), 'svg');
saveas(gcf, fullfile('output', 'consensusBPCsPooled', 'barGraphsConsBPC'), 'png');

% print number of stim sites without BPC
fprintf('%d stim sites with no BPC category\n', sum(isnan(BPCs_expanded(:, 1))));

% print number of stim sites with SNR < 1
fprintf('%d stim sites with SNR < 1\n', sum(BPCs_expanded(:, 2) < 1));

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

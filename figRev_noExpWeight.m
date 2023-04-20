%% Generates Figure Response-1 in response to reviewer #3 in 2nd resubmission, creating BPC outputs without exponential weights
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Main subject-level output file
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

clear; clc;
set(0, 'DefaultFigureRenderer', 'painters');

subNames = {'1', '2', '3', '4', '5'};
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, ...
            {'RB5', 'RB6', 'RC6', 'RC7'}, {'RA9', 'RA10'}};
hemis = 'lrlrr';
subs = struct();
for ii = 1:length(subNames)
    subs(ii).name = subNames{ii};
    subs(ii).chs = subChs{ii};
    subs(ii).hemi = hemis(ii);
end

%% Load mef data for all relevant channels in subject (save as in saveBPCs)

subInd = 1;
sub = subs(subInd).name;
chs = subs(subInd).chs;

subDir = fullfile(pwd, 'data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg');
mkdir(fullfile('output', 'noExp', sprintf('sub-%s', sub)));

mefPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_ieeg.mefd', sub));
channelsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_channels.tsv', sub));
eventsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_events.tsv', sub));

mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'});
mefObj.loadMefTrials([-1, 2]);
mefObj.car(true); % by 64-block
mefObj.pruneChannels(chs);
mefObj.removeLN('SpectrumEstimation');
mefObj.subtractBaseline([-0.5, -0.05], 'median');
mefObj.plotInputs([], [], 200);

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

%% Get data for channel and exclude stim trials at channel

ch = 'LC6'; % hard coded
    
fprintf('Getting data from channel %s\n', ch);

sigdata = squeeze(data(strcmp(mefObj.channels.name, ch), :, :))';
events_ch = events;

% remove stimulated channels from current events
sigdata(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
events_ch(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
sites = groupby(events_ch.electrical_stimulation_site);

% resolve only 4mA or 6mA trials for stim sites with both
for ii = 1:size(sites, 1)
    idxes = sites{ii, 2};
    trialsCurr = events_ch(idxes, :);
    if length(unique(trialsCurr.electrical_stimulation_current)) == 1, continue; end

    trials6ma = idxes(strcmp(trialsCurr.electrical_stimulation_current, '6.0 mA')); % trials at 6 mA stim
    if length(trials6ma) >= 8
        sites{ii, 2} = trials6ma;
        continue
    else % assign the 4.0 mA trials
        sites{ii, 2} = idxes(strcmp(trialsCurr.electrical_stimulation_current, '4.0 mA'));
    end
end
sites(cellfun(@length, sites(:, 2)) < 8, :) = []; % exclude sites with < 8 trials

%% get BPCs

disp('Getting BPCs');
tau = 0.1;
rng('default');

pairTypes = struct('pair', sites(:, 1), 'indices', sites(:, 2));
V = sigdata(:, tt >= 0.011 & tt < 0.5)';
ttBPC = tt(tt >= 0.011 & tt < 0.5);

% NO EXPONENTIAL WEIGHTING

% configure zeta (2 clusters at zeta = 1, 3 clustesr at zeta = 1.35)
zeta = input('zeta threshold = ?\n');
[B, exc, H] = bpc_identify(V, pairTypes, zeta, 50); % 50 reruns

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
plotTrials(ttBPC, curves, 0.2, [], cm, 'LineWidth', 1.5);
saveas(gcf, fullfile('output', 'noExp', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCs', sub, ch)), 'png');
saveas(gcf, fullfile('output', 'noExp', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCs', sub, ch)), 'svg');
close(f)

%% save plotweights

[~, order] = sortrows([BPCs_expanded(:, 1), -BPCs_expanded(:, 2)]); % sort by BPC number and descending plotweight)
T = table(sites(order, 1), BPCs_expanded(order, 1), BPCs_expanded(order, 2), ...
          'VariableNames', {'electrical_stimulation_site', 'BPC', 'SNR'});
writetable(T, fullfile('output', 'noExp', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCSNR.tsv', sub, ch)), ...
          'Filetype','text', 'Delimiter','\t'); 

%% Plot electrodes on gifti brain rendering with BPC SNRs

disp('Plotting Gifti');
gii = gifti(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', sub), sprintf('pial.%s.surf.gii', upper(hemis(subInd)))));
electrodes = readtableRmHyphens(fullfile('data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg', sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub)));
xyzs = [electrodes.x, electrodes.y, electrodes.z];
xyzsPair = ieeg_getPairXyzs(split(sites(:, 1), '-', 2), electrodes);

f = figure('Position', [1000, 100, 1000, 800]); % normal gifti

ieeg_RenderGifti(gii); alpha 0.2; hold on
switch hemis(subInd)
    case 'r'
        ieeg_viewLight(90, -40);
    case 'l'
        ieeg_viewLight(-90, -40);
end

plot3(xyzs(:,1), xyzs(:,2), xyzs(:,3), 'o', 'Color', 'k', 'MarkerSize', 6, 'MarkerFaceColor', 'w');
%text(xyzs(:,1), xyzs(:,2), xyzs(:,3), electrodes.name, 'Color', 'k');
tgt = find(strcmp(electrodes.name, ch)); % circle target electrode
plot3(xyzs(tgt,1), xyzs(tgt,2), xyzs(tgt,3), 'o', 'Color', 'r', 'MarkerSize', 12, 'LineWidth', 3);

for b = 1:max(BPCs_expanded(:, 1))
    ix_bool = BPCs_expanded(:, 1) == b;
    color_add_custom(xyzsPair(ix_bool, :), BPCs_expanded(ix_bool, 2), cm(b, :), 0.8, 3, [6, 20], 's');
end
color_add_custom(xyzsPair(isnan(BPCs_expanded(:, 1)), :), 0.5, 0.1*[1 1 1], 0.8, 1, [6, 10], 's'); %non-BPCs
hold off

saveas(f, fullfile('output', 'noExp', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_giftiBPCs', sub, ch)), 'png');
set(f, 'Position', [1000, 100, 600, 500])
saveas(f, fullfile('output', 'noExp', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_giftiBPCs_small', sub, ch)), 'png');
close(f)

%% Plot mean trial of each stim site for each BPC

disp('Plotting mean traces of stim sites');
meanTrialDir = fullfile('output', 'noExp', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_meanTrials', sub, ch));
mkdir(meanTrialDir);

for kk=1:length(B)
    f = plot_meanTrials(tt, sigdata, sites(B(kk).pairs, :), 200, 'LineWidth', 1, 'Color', cm(kk, :));
    xlim([0, 0.5]);
    xline(0.011, 'Color', 'r');

    labs = flip(get(gca, 'YTickLabel'));
    labs(B(kk).plotweights > 1) = cellfun(@(s) sprintf('*%s', s), labs(B(kk).plotweights > 1), 'UniformOutput', false);
    set(gca, 'YTickLabel', flip(labs));

    saveas(f, fullfile(meanTrialDir, sprintf('sub-%s_ch-%s_BPC%d_meanTrials', sub, ch, kk)), 'png');
    saveas(f, fullfile(meanTrialDir, sprintf('sub-%s_ch-%s_BPC%d_meanTrials', sub, ch, kk)), 'svg');
    close(f)
end

if ~isempty(exc)
    f = plot_meanTrials(tt, sigdata, sites(exc, :), 200, 'LineWidth', 1, 'Color', 'k');
    xlim([0, 0.5]);
    xline(0.011, 'Color', 'r');

    saveas(f, fullfile(meanTrialDir, sprintf('sub-%s_ch-%s_exc_meanTrials', sub, ch)), 'png');
    saveas(f, fullfile(meanTrialDir, sprintf('sub-%s_ch-%s_exc_meanTrials', sub, ch)), 'svg');
    close(f)
end



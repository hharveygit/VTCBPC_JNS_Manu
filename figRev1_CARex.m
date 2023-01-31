%% Code to plot example of CAR channels used to rereference a CCEP trial
%   2022/12/19
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
%% Configure subject and load CCEP data

clear;
set(0, 'DefaultFigureRenderer', 'painters');

sub = '1';
chs = 'LC6';

subDir = fullfile(pwd, 'data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg');

mefPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_ieeg.mefd', sub));
channelsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_channels.tsv', sub));
eventsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_events.tsv', sub));

mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'});
mefObj.loadMefTrials([-1, 2]);

% extract data
tt = mefObj.tt;
srate = mefObj.srate;
events = mefObj.evts;
channels = mefObj.channels;
data = mefObj.data;

%% Choose CCEP trial, get line noise removed version

tr = 74; % trial index
dataTr = data(:, :, tr);

% remove line noise and baseline for chosen trial
dataTrClean = removeLineNoise_SpectrumEstimation(dataTr, srate, 'LF = 60, NH = 3, HW = 3', false);
dataTrClean = dataTrClean - median(dataTrClean(:, tt >= -0.5 & tt < -0.05), 2);

%% Determine included and excluded CAR channels

seegChs = strcmp(channels.type, 'SEEG'); % all SEEG channels

badChs = ~strcmp(channels.status, 'good') & seegChs; % not marked good
stimChs = ismember(channels.name, split(events.electrical_stimulation_site(tr), '-')); % stim channels

% metrics to filter on
varLate = var(dataTr(:, tt > 0.5 & tt < 1), [], 2);
varEarly = var(dataTr(:, tt > 0.015 & tt < 0.1), [], 2);

varLateThresh = prctile(varLate(seegChs & ~badChs & ~stimChs), 95); % 95% noise threshold for non-stim and non-bad channels
minEarlyVar = prctile(varEarly(seegChs & ~badChs & ~stimChs), 75); % assume bottom 75%tile of channels don't have early responses

% these two may overlap
noisyChs = varLate > varLateThresh & seegChs & ~badChs & ~stimChs;
respChs = varEarly ./ varLate > 1.5 & varEarly > minEarlyVar & seegChs & ~badChs & ~stimChs;

goodChs = seegChs & ~badChs & ~stimChs & ~noisyChs & ~respChs; % good channels to be included in CAR

%% Make plot of all channels and highlight different categories.

mkdir(fullfile('output', 'CARsensitivityTest'));

% plot good channels in 2 separate sets for ease of visualization
idxGC = find(goodChs);
set1 = idxGC(1:ceil(length(idxGC)/2));
set2 = idxGC(ceil(length(idxGC)/2)+1:end);

figure('Position', [200, 200, 400, length(set1)*15]); plotTrials(tt, dataTr(set1, :)', 500, set1, 'k');
xlim([0, 1]); title('Good Chs, set 1'); xline(0.015, 'Color', 'r'); xline(0.1, 'Color', 'r');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_goodChs1', sub, chs, tr)), 'png');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_goodChs1', sub, chs, tr)), 'svg');

figure('Position', [200, 200, 400, length(set2)*15]); plotTrials(tt, dataTr(set2, :)', 500, set2, 'k');
xlim([0, 1]); title('Good Chs, set 2'); xline(0.015, 'Color', 'r'); xline(0.1, 'Color', 'r');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_goodChs2', sub, chs, tr)), 'png');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_goodChs2', sub, chs, tr)), 'svg');

% Plot stimulated channels
figure('Position', [200, 200, 400, 200]); plotTrials(tt, dataTr(stimChs, :)', 2500, find(stimChs), 'k');
title('Stimulated channels (excluded)');
xlim([0, 1]); %title('Stimulated channels');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_stimChs', sub, chs, tr)), 'png');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_stimChs', sub, chs, tr)), 'svg');

% Plot noisy and responsive channels together
noisyNresp = [find(noisyChs & ~respChs); find(noisyChs & respChs); find(respChs & ~noisyChs)]; % order to make highlighting easier
figure('Position', [200, 200, 400, length(noisyNresp)*20]);
ys = plotTrials(tt, dataTr(noisyNresp, :)', 500, noisyNresp, 'k');
title('Noisy and responsive channels (excluded)');
xlim([0, 1]); %title('Noisy and responsive channels');
xline(0.015, 'Color', 'r'); xline(0.1, 'Color', 'r');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_noisyNrespChs', sub, chs, tr)), 'png');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_noisyNrespChs', sub, chs, tr)), 'svg');

%% Save plots of the actual common average references (corresponding to each headbox)

cars = [];
for ii = 1:4
    toUse = idxGC(idxGC > (ii-1)*64 & idxGC <= ii*64);
    if isempty(toUse), break; end
    
    cars = [cars; mean(dataTr(toUse, :))];
end

% each line is the CAR for a 64-channel headbox
figure('Position', [200, 200, 400, 300]);
plotTrials(tt, cars', 500, [], 'k');
xlim([0, 1]); title('CARs');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_CARs', sub, chs, tr)), 'png');
saveas(gcf, fullfile('output', 'CARsensitivityTest', sprintf('sub-%s_%s_tr=%d_CARs', sub, chs, tr)), 'svg');

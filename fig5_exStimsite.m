%% Used to plot the evoked potentials in Figure 5A (bottom).
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Creates example CCEP plot for figure 5
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
%% Load data for subject 2

clear;
set(0, 'DefaultFigureRenderer', 'painters');

sub = '2';
ch = 'RC6';
subDir = fullfile(pwd, 'data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg');

mefPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_ieeg.mefd', sub));
channelsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_channels.tsv', sub));
eventsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_events.tsv', sub));

mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'});
mefObj.loadMefTrials([-1, 2]);
mefObj.car(true); % by 64-block
mefObj.pruneChannels(ch);
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

% remove stimulated channels from current events
sigdata = squeeze(data)';
events_ch = events;
sigdata(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
events_ch(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
sites = groupby(events_ch.electrical_stimulation_site);

%% Plot stim site and save

% Remove stim artifact with forward/back artifact removal method, before wavelet or broadband filtering
sigdataNoart = remove_artifact(sigdata', find(tt >= -11e-3 & tt <= 11e-3), 'crowther')';

stimsite = 'RC2-RC3';
sigdataStim = sigdataNoart(sites{strcmp(sites(:, 1), stimsite), 2}, :);

figure('Position', [200, 200, 600, 300]); hold on
plot(tt, sigdataStim', 'Color', [0.3, 0.3, 0.3]);
plot(tt, mean(sigdataStim), 'k', 'LineWidth', 1.5);
yline(0, 'Color', [0.5, 0.5, 0.5]);
xlim([-0.1, 1]); ylim([-2000, 1000]);
saveas(gcf, fullfile('output', 'spectralBB', 'fig5AExCcep'), 'png');
saveas(gcf, fullfile('output', 'spectralBB', 'fig5AExCcep'), 'svg');

%% Code to compare line noise removal methods
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
%% Load subject data

clear;
set(0, 'DefaultFigureRenderer', 'painters');

sub = '1';
chs = {'LC6', 'LB6', 'LB7', 'LC5'};

subDir = fullfile(pwd, 'data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg');

mefPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_ieeg.mefd', sub));
channelsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_channels.tsv', sub));
eventsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_events.tsv', sub));

mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'});
mefObj.loadMefTrials([-1, 2]);
mefObj.car(true); % by 64-block
mefObj.pruneChannels(chs);
mefObj.subtractBaseline([-0.5, -0.05], 'median');

dataNoisy = mefObj.data; % exctract data without any line noise removal

% copy of mef object to be processed with notch filter
mefObjNotch = copy(mefObj);

mefObj.removeLN('SpectrumEstimation');
mefObj.subtractBaseline([-0.5, -0.05], 'median');

mefObjNotch.removeLN('notch');
mefObjNotch.subtractBaseline([-0.5, -0.05], 'median');

% extract data
tt = mefObj.tt;
srate = mefObj.srate;
events = mefObj.evts;
channels = mefObj.channels;

dataInterpol = mefObj.data; % data with line noise removed by spectrum interpolation
dataNotch = mefObjNotch.data; % data with line noise removed by notch filter

%% Get data for one stim site

sites = groupby(events.electrical_stimulation_site);

% choose channel and stim site
ch = 1;
stimIdx = 12;

dataNoisyTr = squeeze(dataNoisy(ch, :, sites{stimIdx, 2}));
dataInterpolTr = squeeze(dataInterpol(ch, :, sites{stimIdx, 2}));
dataNotchTr = squeeze(dataNotch(ch, :, sites{stimIdx, 2}));

mkdir(fullfile('output', 'LNremoval'));

% Plot time series examples
% no line noise removed
figure('Position', [200, 200, 800, 400]); plot(tt, dataNoisyTr, 'Color', 'k'); yline(0, 'Color', [0.5, 0.5, 0.5]); xlim([-0.1, 0.5]); ylim([-200, 600]);
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoiseExNoisy_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'png');
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoiseExNoisy_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'svg');

% line noise removed with spectrum interpolation (used in this manuscript)
figure('Position', [200, 200, 800, 400]); plot(tt, dataInterpolTr, 'Color', 'k'); yline(0, 'Color', [0.5, 0.5, 0.5]); xlim([-0.1, 0.5]); ylim([-200, 600]);
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoiseExInterpol_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'png');
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoiseExInterpol_sub%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'svg');

% line noise removed with notch filter
figure('Position', [200, 200, 800, 400]); plot(tt, dataNotchTr, 'Color', 'k'); yline(0, 'Color', [0.5, 0.5, 0.5]); xlim([-0.1, 0.5]); ylim([-200, 600]);
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoiseExNotch_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'png');
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoiseExNotch_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'svg');


% Plot welch's estimate PSDs
[pxxNoisyTr, f] = pwelch(dataNoisyTr(tt >= 0.011 & tt <= 0.5, :), hann(srate/8), srate/16, 0:200, srate);
pxxInterpolTr = pwelch(dataInterpolTr(tt >= 0.011 & tt <= 0.5, :), hann(srate/8), srate/16, 0:200, srate);
pxxNotchTr = pwelch(dataNotchTr(tt >= 0.011 & tt <= 0.5, :), hann(srate/8), srate/16, 0:200, srate);

% no line noise removed
figure('Position', [200, 200, 800, 400]); plot(f, log10(pxxNoisyTr), 'Color', 'k'); ylim([-4, 4]);
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoisePxxNoisy_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'png');
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoisePxxNoisy_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'svg');

% line noise removed with spectrum interpolation
figure('Position', [200, 200, 800, 400]); plot(f, log10(pxxInterpolTr), 'Color', 'k'); ylim([-4, 4]);
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoisePxxInterpol_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'png');
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoisePxxInterpol_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'svg');

% line noise removed with notch filter
figure('Position', [200, 200, 800, 400]); plot(f, log10(pxxNotchTr), 'Color', 'k'); ylim([-4, 4]);
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoisePxxNotch_sub-%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'png');
saveas(gcf, fullfile('output', 'LNremoval', sprintf('lineNoisePxxNotch_sub%s_%s_%s', sub, chs{ch}, sites{stimIdx, 1})), 'svg');


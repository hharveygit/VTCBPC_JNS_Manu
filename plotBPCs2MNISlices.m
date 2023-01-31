%% This script is used after clusterBPCs.m to plot stim sites to MNI 152 T1 MRI slices
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Plots stim sites to MNI152 T1 MRI slices, colored by consensus BPC
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
%% Configure paths to MNI 152 output, MNI 152 nifti

clear;
set(0, 'DefaultFigureRenderer', 'opengl');

mkdir(fullfile('output', 'clusterBPCs', 'MNI152slices'));

niiMNI = fullfile('data', 'derivatives', 'mni152.nii'); % MRI152 nifti from MRIcroGL
masterTab = readtable(fullfile('output', 'clusterBPCs', 'masterMNItab.tsv'), 'FileType', 'text', 'Delimiter', '\t');

% Subs 2 and 4 are in right hemisphere. Coordinates were inverted to left when stored to masterMNItab.tsv. Here we invert back for plotting
originalX = masterTab.x;
originalX(ismember(masterTab.sub, [2, 4])) = -originalX(ismember(masterTab.sub, [2, 4]));

cmap = getCmapVTC('bpc');
cmap(5, :) = [0.3, 0.3, 0.3];
         
%% Load interactive GUI with MNI T1

electrodes = table({'dummy'}, 0, 0, -50, 'VariableNames', {'name', 'x', 'y', 'z'}); % fake electrodes input necessary for sliceGUI to run

% voxel dim = 0.7375 mm
% make cuts every 8.112 mm (11 voxels)
opts.plotnames = false;
opts.elecWidth = 8.112/2; % half of the cut jump is used as width to snap electrodes too
opts.initialPos = 0;

hands = sliceGUI(niiMNI, electrodes, opts); % open up the slice plot GUI

%% Calculate range of slices to make

% most inferior and superior electrode positions across subjects 1-4 in MNI152 coordinates 
disp([min(masterTab.z), max(masterTab.z)]);

% make cuts every 8.112 mm (11 voxels) to span all electrodes in the axial dimension
% These are the z-positions used in Figure 4.
coronalSliceCenters = [-32.7, -24.6, -16.5, -8.3, -0.2, 7.9, 16.0, 24.1, 32.2, 40.3, 48.4, 56.6];

% In the slice plot GUI, click on the coronal slice at each z-level (shown in bottom left corner of axial slice),
% before running the section below to plot stim sites at that z-level
% Then click on a different z-level, and run section below again, repeat until all desired slices have been saved.

% pause here if code is run automatically
return;

%% Add stim sites across subjects 1-4, colored by consensus BPC
% To plot a different slice, click somewhere else in the GUI and run this section again

locs = [originalX, masterTab.y, masterTab.z];
wts = [masterTab.SNR, masterTab.consBPC];

optsNan.marker = 's';
optsNan.szrng = [4, 8];
addSliceStims(hands, locs(isnan(masterTab.SNR) | masterTab.SNR < 1, :), [], [0.1, 0.1, 0.1], optsNan); % non-assigned sites or SNR < 1

optsStims.clip = 3;
optsStims.marker = 's';
optsStims.szrng = [8, 15];
optsStims.cscale = 0.6;
addSliceStims(hands, locs(masterTab.SNR >= 1, :), wts(masterTab.SNR >= 1, :), cmap, optsStims);

% save ong and eps (kjm_printfig locate in mnl_ieegBasics/functions)
kjm_printfig(fullfile('output', 'clusterBPCs', 'MNI152slices', sprintf('MNI152BPCs_%s', datestr(now, 'yymmdd-HHMMSS'))), [30, 10]);

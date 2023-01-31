%% Plots stim sites to subject T1 MRI slices, colored by consensus BPC
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Plots stim sites to subject T1 MRI slices, colored by consensus BPC
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
%% Configure subjects, paths

clc, clear;
set(0, 'DefaultFigureRenderer', 'opengl');

subNames = {'1', '2', '3', '4'};
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6'}};
        
% Load the ACPC-space nifti path and corresponding ACPC electrode positions for each subject
niiAcpc = cell(4, 1);
elecPaths = cell(4, 1);
for ii = 1:4
    switch ii
        case 1 % sub1 T1w MRI and electrodes used in the rest of this analysis are already in AC-PC space
            niiAcpc{ii} = fullfile('data', 'derivatives', 'freesurfer', 'sub-1', 'sub-1_ses-mri01_T1w_acpc_deFaced.nii');
            elecPaths{ii} = fullfile('data', 'sub-1', 'ses-ieeg01', 'ieeg', 'sub-1_ses-ieeg01_electrodes.tsv');
        otherwise
            niiAcpc{ii} = fullfile('data', 'derivatives', 'acpc', sprintf('sub-%d', ii), sprintf('sub-%d_ses-mri01_T1w_acpc_deFaced.nii', ii));
            elecPaths{ii} = fullfile('data', sprintf('sub-%d', ii), 'ses-ieeg01', 'ieeg', sprintf('sub-%d_ses-ieeg01_space-ACPC_electrodes.tsv', ii));
    end
end

cmap = getCmapVTC('bpc');
cmap(5, :) = [0.3, 0.3, 0.3]; % dark gray for excluded BPC stim sites
         
%% Load interactive GUI for specific subject. Repeat this section for each subject

sub = input('\nSubject? Choose from 1 - 4\n');

electrodes = readtableRmHyphens(elecPaths{sub});

opts.plotnames = false;
opts.elecWidth = 8;
opts.initialPos = 0;

if ismember(sub, [2, 4])
    opts.clim = [0, 0.6];
elseif sub == 3
    opts.clim = [0, 1];
else
    opts.clim = [0, 0.8];
end

hands = sliceGUI(niiAcpc{sub}, electrodes, opts);

return % pause here to read instructions below

% Manually click on the coronal slice at the (Z, X) coordinates desired for the current subject before running the next section.
% (Z, X) coordinates plotted in Figure 7 are as follows for each subject:
% S1 = (-15.0, -23.0), S2 = (-15.0, 23.0), S3 = (-13.0, -24.0), S4 = (-13.0, 26.5)

% To repeat for multiple measurement electrodes in the same subject, click anywhere else in the coronal slice to clear current stim sites plotted

%% Add stimulation sites for the desired measurement electrode

ch = input(sprintf('Electrode number? Choose from 1 - %d for current subject\n', length(subChs{sub})));
mkdir(fullfile('output', sprintf('sub-%s', subNames{sub}), 'hippAmygSlices'));

bpcTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', subNames{sub}), sprintf('sub-%s_ch-%s_BPCSNR.tsv', subNames{sub}, subChs{sub}{ch})), 'electrical_stimulation_site', 1);
pairLocs = ieeg_getPairXyzs(split(bpcTable.electrical_stimulation_site, '-', 2), electrodes);

% non-assigned stim sites
optsNan.marker = 's';
optsNan.szrng = [3, 7];
addSliceStims(hands, pairLocs(isnan(bpcTable.SNR), :), [], [0.1, 0.1, 0.1], optsNan)

optsStims.clip = 3;
optsStims.marker = 's';
optsStims.szrng = [6, 12];
optsStims.cscale = 0.6;
addSliceStims(hands, pairLocs, [bpcTable.SNR, bpcTable.consBPC], cmap, optsStims);

% Save
kjm_printfig(fullfile('output', sprintf('sub-%s', subNames{sub}), 'hippAmygSlices', sprintf('sub-%s_ch-%s_slice', subNames{sub}, subChs{sub}{ch})), [30, 10]);
fprintf('Stim sites plotted and saved to T1 MRI slice for subject %s, ch %s\n', subNames{sub}, subChs{sub}{ch});

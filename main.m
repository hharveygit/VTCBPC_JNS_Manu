%% Use each section of this main script to output all results and figures for the manuscript.
%   2022/01/31
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: main script.
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
%% Path configuration

% Make sure you are in the main manuscript directory (VTCBPC)
% Data should be downloaded from Openneuro and copied to a folder here (./data)

% Add paths to dependencies
addpath(genpath('path/to/mnl_ieegBasics'));
addpath(genpath('path/to/matmef'));
addpath('path/to/mnl_seegview');
addpath(genpath('path/to/vistasoft'));
addpath('path/to/spm12'); addpath(fullfile('path/to/spm12', 'config'));
addpath(genpath('path/to/freesurfer/7.1.1/matlab'));
addpath(genpath('path/to/gifti'));

addpath(genpath('functions'));

%% 1a. Calculate MNI 152 forward deformation transformations using SPM 12

spm fmri

% Press "Batch"
% SPM > Spatial > Segment
% Specify the subject's T1w nifti file for the "Volumes" input.
% Select "Forward" for "Deformation Fields".
% Press the play button
% Keep the file titled "y_[t1FileName].nii". all other output files can be deleted

%% 1b. Use forward deformation transformations to calculate MNI152 electrode positions for subjects 1-4

% Output is saved to 'output\sub-X\MNI152\.' for each subject X
saveMNIelectrodes

%% 2. Calculate number of stim sites for each subject - table 1 values

determineStimSites

%% 3. Creates methods figure 1 panels C, D, and E

% Output is saved to 'output\methodsFigure\.'
methodfigMatrices

%% 4. For figures 2 and 6. Save subject-level BPCs and per-subject wavelet spectrogram, broadband estimates, for all measurement elecs.

% Output is saved to 'output\sub-X\.' for each subject X
saveBPCs

%% 5. Save Destrieux cortical labels to each stim site in the output BPCSNR.tsv files.

assignStimDestrieux
% LA1-LA2 stim site in sub-1 was manually corrected in tsv output to "Left_Amygdala" by visual inspection of its position in CT/T1w MRI overlay.

%% 6. Calculate consensus BPCs and anatomical distribution of stim sites color coded by consensus BPCs on brain, plus cross-val

% Output is saved to 'output\clusterBPCs\.'
clusterBPCs

% Save legend of cortical and subcortical structures (brains in Figure 3D, right), to 'output\anatGifti_sub3_cortical_*.png'
bpcRegions

%% 7. Plot stim sites colored by consensus BPCs on slices of the MNI152 T1w MRI

% You will need to manually click at different axial z-levels of the GUI, and then running the last section of the script to save
% stim sites plotted at that z-level.
% Output is saved to 'output\clusterBPCs\MNI152slices\.'
plotBPCs2MNISlices

%% 8. Calculate spectrograms and broadband estimates for consensus BPCs, across measurement electrodes

% Output is saved to 'output\spectralBBs\.'
bpcSpectraBB

% Save the subject 2 example CCEP in figure 5A (top)
% Output is saved to 'output\spectralBBs\.'
fig5_exStimsite

%% 9. Plot stim sites by consensus BPC on inflated subject pial renderings. Also plot BPCs colored/ordered by consensus identity

% Output is saved to 'output\sub-X\inflated\.'
plotInflatedBPCs

%% 10. Plot stim sites by consensus BPC color on subject T1 MRI slices

% Sections 2 and 3 of plotBPCs2Slice must be run iteratively for each subject and measurement electrode number:
% Running section 2 pulls up the interactive T1 MRI slice GUI. Click on the (Z, X) positions in the coronal slice
% corresponding to each subject before running section 3.
% (Z, X) coordinates plotted in Figurei 7 are as follows for each subject:
% S1 = (-15.0, -23.0), S2 = (-15.0, 23.0), S3 = (-13.0, -24.0), S4 = (-13.0, 26.5)

% Outputs are saved to 'output\sub-X\hippAmygSlices\.'
plotBPCs2Slice

%% 11. Create Basis Profile Spectrograms (BPS) for all subjects/measurement electrodes, and quantify similarity to BPCs

% Outputs are saved to 'output\BPS\.'
overlapWithBPS

%% 12. Monopolar subject-level BPC outputs without common average rereferencing (CAR sensitivity analysis)

% Outputs are saved to 'output\CARsensitivityTest\sub-1\.'
BPCsNoCAR

% Run this to generate outputs of channels that are included/excluded in the generation of the CAR for a CCEP trial
% Outputs are saved to 'output\CARsensitivityTest\.'
figRev1_CARex

%% 13. Comparison of the spectrum interpolation line noise removal used in this manuscript to standard notch filter

% Outputs are saved to 'output\LNremoval\.'
figRev2_LNremoval

%% 14. Distance between stim site - measurement electrode vs. consensus BPC assignment or BPC SNR

% Outputs are saved to 'output\distanceVBPC\.'
distanceVBPC

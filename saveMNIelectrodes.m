%% saveMNIelectrodes.m. Saves MNI 152 coordinates for subject electrodes using SPM12.
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Save MNI 152 electrode coordinates using SPM12 outputs
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
%% Iterate across subjects

clear;

for sub = 1:4
    %% Configure paths to electrodes, T1w MRI, forward deformation files
    fprintf('sub-%d:\n', sub);
    
    % Path to electrodes file
    elecPath = fullfile('data', sprintf('sub-%d', sub), 'ses-ieeg01', 'ieeg', sprintf('sub-%d_ses-ieeg01_electrodes.tsv', sub));
    
    % Paths to T1w MRI and forward deformation files (generated from SPM12). Subject 1 uses AC-PC MRI; subjects 2-4 use native space
    if sub == 1
        t1Path = fullfile('data', 'derivatives', 'freesurfer', 'sub-1', 'sub-1_ses-mri01_T1w_acpc_deFaced.nii');
        defPath = fullfile('data', 'derivatives', 'freesurfer', 'sub-1', 'y_sub-1_ses-mri01_T1w_acpc.nii');
    else
        t1Path = fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%d', sub), sprintf('sub-%d_ses-mri01_T1w_deFaced.nii', sub));
        defPath = fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%d', sub), sprintf('y_sub-%d_ses-mri01_T1w.nii', sub));
    end
    
    % Location to save the electrode images generated + MNI 152 electrodes file
    mniOutDir = fullfile('output', sprintf('sub-%d', sub), 'MNI152');
    mkdir(fullfile('output', sprintf('sub-%d', sub), 'MNI152'));
    
    %% Read electrodes and use SPM12 to convert to MNI152 space, save.
    
    electrodes = readtableRmHyphens(elecPath);

    xyz = [electrodes.x, electrodes.y, electrodes.z];
    xyzMni = ieeg_getXyzMniManu(xyz, t1Path, defPath, mniOutDir);

    electrodesMni = electrodes;
    electrodesMni.x = xyzMni(:, 1);
    electrodesMni.y = xyzMni(:, 2);
    electrodesMni.z = xyzMni(:, 3);
    
    writetable(electrodesMni, fullfile(mniOutDir, sprintf('sub-%d_ses-ieeg01_space-MNI152NLin6Sym_electrodes.tsv', sub)), 'Filetype', 'text', 'Delimiter', '\t'); 
    
end
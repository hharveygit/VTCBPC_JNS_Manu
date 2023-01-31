%% Save column of Destrieux labels to BPCSnr tsv files output by saveBPCs
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: assigns anatomical labels to stim sites in BPC SNR outputs
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

% directory to subject-level Freesurfer outputs
clear;
FSDir = fullfile('data', 'derivatives', 'freesurfer');
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6', 'RC6', 'RC7'}, {'RA9', 'RA10'}};
        
% iterate across the 5 subjects
for ii = 1:5
    
    sub = num2str(ii);
    for jj = 1:length(subChs{ii})
        
        ch = subChs{ii}{jj};
        fprintf('sub-%s %s\n', sub, ch);
        
        % load electrode positions
        electrodes = readtableRmHyphens(fullfile('data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg', sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub)));
        
        % load tsv output from saveBPCs
        snrTablePath = fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCSNR.tsv', sub, ch));
        snrTable = readtableRmHyphens(snrTablePath, 'electrical_stimulation_site', 1);
        
        stimPairs = split(snrTable.electrical_stimulation_site, '-', 2);
        pairXyzs = ieeg_getPairXyzs(stimPairs, electrodes); 
        labelsDestrieux = ieeg_getLabelXyzDestrieux(pairXyzs, fullfile(FSDir, sprintf('sub-%s', sub)));

        snrTable.labelDestrieux = labelsDestrieux;
        writetable(snrTable, snrTablePath, 'Filetype','text', 'Delimiter','\t'); % write back to same location
    end
    
end
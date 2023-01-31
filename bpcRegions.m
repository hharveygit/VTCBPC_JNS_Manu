%% Create cortical and subcortical legends in figure 3 (consensus BPC)
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%    VTCBPC manuscript package: Creates anatomical legends for cortical and
%    subcortical regions
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
%% Render gifti with areas colored by anatomical colormap

clear;
set(0, 'DefaultFigureRenderer', 'painters');

% Use subject 3 as example brain for the cortical and subcortical labels
annotPath = fullfile('data', 'derivatives', 'freesurfer', 'sub-3', 'label', 'lh.aparc.a2009s.annot'); % anatomical segmentation
gii = gifti(fullfile('data', 'derivatives', 'freesurfer', 'sub-3', 'pial.L.surf.gii'));
                
cmAnat = getCmapVTC('anat'); % colormap of different anatomical regions
[~, label, colortable] = read_annotation(annotPath);
label = changem(label, 1:size(colortable.table(:,5)), colortable.table(:,5)); % change labels to be 1:n index
label_cmap = 0.8*ones(size(colortable.table, 1), 3); % initialize to gray; matches cmAnat(6, :)

anatBool = @(labs) contains(lower(colortable.struct_names), labs); % returns boolean of where labels are found in colortable
label_cmap(anatBool({'oc-temp', 'collat'}), :) = repmat(cmAnat(3, :), sum(anatBool({'oc-temp', 'collat'})), 1);
label_cmap(anatBool({'temporal_', 'temp_sup'}), :) = repmat(cmAnat(4, :), sum(anatBool({'temporal_', 'temp_sup'})), 1);
label_cmap(anatBool('ins'), :) = repmat(cmAnat(5, :), sum(anatBool('ins')), 1);

f = figure('Position', [800, 600, 800, 600]);
ieeg_RenderGiftiLabels(gii, label, label_cmap, colortable.struct_names);
ieeg_viewLight(-90*1.1, -30);
saveas(f, fullfile('output', 'anatGifti_sub3_cortical'), 'png');


%% Render subcortical structures (amygdala & hippocampus)

% subcortical segmentation
[segs, M] = load_mgh(fullfile('data', 'derivatives', 'freesurfer', 'sub-3', 'mri', 'aseg.mgz'));

% labels corresponding to segmentation structures to render; cortex = 3, left hippocampus = 17, left amygdala = 18
numLabs = [17, 18];
cmAnat = getCmapVTC('anat');

f = figure('Position', [800, 600, 1000, 800]); hold on
ieeg_RenderGifti(gii); alpha 0.3
renderSubcortical(numLabs, segs, M, cmAnat);
hold off
ieeg_viewLight(-90*1.1, -30);
saveas(f, fullfile('output', 'anatGifti_sub3_subcortical'), 'png');


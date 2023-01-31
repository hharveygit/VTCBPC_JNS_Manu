%% Renders subcortical structures from freesurfer segmentation
%
%   renderSubcortical(labels, segs, M);
%   renderSubcortical(labels, segs, M, cmap);
%       labels =        1xa double, numerical label (or list of labels) of structures to render. Corresponds to labels
%                           in mri/talairach.label_intensities.txt.
%                           E.g. L hippocampus = 17, L amygdala = 18, R hippocampus = 53, R amygdala = 54
%       segs =          nxmxp double, segmentation volume loaded from load_mgh. Values of volume correspond to identity
%                           of voxels
%       M =             4x4 double, transformation matrix to go from voxel -> coordinates, loaded from load_mgh
%       cmap =          bx3 double [0...1], (optional). Colormap to use for structures, where each row is a color. 1st
%                           row = color of 1st structure in <labels>, 2nd row = color of 2nd structure, etc...
%                           Default Matlab colormap if not given.
%
%   Example usage:
%
%       % Make sure freesurfer functions are added to path
%       [segs, M] = load_mgh('mri/aseg.mgz');       % load segmentations
%       labels = [17, 18];                          % left hippocampus & amygdala
%       figure; hold on
%       ieeg_RenderGifti(gifti); alpha 0.3          % render transparent cortex
%       renderSubcortical(labels, segs, M);         % render hipp, amyg, using default matlab colormap
%       hold off
%
%   HH 2021
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
function renderSubcortical(labels, segs, M, cmap)

    if nargin < 4
        cmap = get(0,'defaultAxesColorOrder'); % default colormap
    end

    for ii = 1:length(labels)

        idx = [];
        [idx(:, 1), idx(:, 2), idx(:, 3)] = ind2sub(size(segs), find(segs == labels(ii))); % 3D voxel indices of structure
        p = [idx - 1, ones(size(idx, 1), 1)]*M'; p(:, 4) = []; % transform to position coordinates

        t = boundary(p); % boundary of point cloud

        trisurf(t, p(:,1), p(:,2), p(:,3), 'facecolor', cmap(ii, :), 'EdgeColor', 'none'); material dull

    end
    
end
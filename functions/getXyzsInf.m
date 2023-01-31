%% This function returns inflated gifti XYZ positions from an electrodes.tsv table
% Adapted from DH ieeg_snap2inflate.m to perform separately for each hemisphere
%
%   xyzsInf = getXyzsInf(electrodes, hemi, gii, giiInf);
%   xyzsInf = getXyzsInf(electrodes, hemi, gii, giiInf, distThresh);
%       electrodes =        nx_ table, read from an electrodes.tsv file. Must contain columns "hemisphere" (str, 'R' or 'L')
%                               and "Destrieux_label" (numerical, where 0 = white matter).
%       hemi =              char, 'r' or 'l', corresponding to the hemisphere of the pial and inflated giftis
%       gii =               gifti object, pial surface of one hemisphere
%       giiInf =            gifti object, inflated surface of one hemisphere
%       distThresh =        double (optional), Only electrodes within <distThresh> mm of a pial vertex are kept in output.
%                               Default = 5.
%
%   Returns:
%       xyzsInf =           nx3 double, XYZ positions of gray matter electrodes in inflated gifti space. Electrodes in
%                               white matter, contralateral hemisphere, or beyond <distThresh> are returned as [NaN, NaN, NaN] rows.
%
%   HH 2021
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Calculates inflated positions from non-inflated pial surface
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
function xyzsInf = getXyzsInf(electrodes, hemi, gii, giiInf, distThresh, checks)

    if nargin < 6, checks = true; end % to perform checks?
    if nargin < 5 || isempty(distThresh), distThresh = 5; end % maximum distance allowed between valid electrode and gifti vertex
    
    hemi = lower(hemi);
    assert(ismember(hemi, {'r', 'l'}), "hemi must be given as 'r' or 'l'");
    assert(size(gii.vertices, 1) == size(giiInf.vertices, 1), 'gii and giiInf do not have the same number of vertices');
    
    xyzs = [electrodes.x, electrodes.y, electrodes.z];
    hemiLabs = lower(electrodes.hemisphere);
    anatLabs = electrodes.Destrieux_label; % numerical labels
    anatText = electrodes.Destrieux_label_text; % char labels

    xyzsInf = nan(size(xyzs));
    for ii = 1:size(xyzs, 1)
        
        if ~strcmp(hemiLabs(ii), hemi) && checks % wrong hemisphere
            continue
        end
        
        if (isnan(anatLabs(ii)) || anatLabs(ii) == 0 || anatLabs(ii) == 41) && checks % non existent, outside brain, or WM
            continue
        end
        
        if ~startsWith(anatText{ii}, sprintf('%sh', hemi)) && checks % must start with lh or rh to be a cortical site
            continue
        end
        
        dists = vecnorm(gii.vertices - xyzs(ii, :), 2, 2); % distance from each vertex in gifti to current electrode
        [minDist, idx] = min(dists);
        
        if minDist <= distThresh
            xyzsInf(ii, :) = giiInf.vertices(idx, :);
        end
        
    end
    
end
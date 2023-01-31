%% This function returns inflated gifti XYZ positions for input gifti XYZ positions
% Adapted from getXyzsInf.m; does not require electrodes.tsv input, so checks positions based on individually-input destrieuxLabels
%
%   xyzsInf = getXyzsInfSimple(xyzs, hemi, gii, giiInf);
%   xyzsInf = getXyzsInfSimple(xyzs, hemi, gii, giiInf, distThresh, destrieuxLabels);
%       xyzs =              nx3 matrix of positions. Rows are [x, y, z] positions to be transformed
%       hemi =              char, 'r' or 'l', corresponding to the hemisphere of the pial and inflated giftis
%       gii =               gifti object, pial surface of one hemisphere
%       giiInf =            gifti object, inflated surface of one hemisphere
%       distThresh =        double (optional), Only electrodes within <distThresh> mm of a pial vertex are kept in output.
%                               Default = 3.
%       destrieuxLabels =   nx1 cell array, corresponding to each position in xyzs. If given, positions in the wrong hemisphere or "Cerebral_White_Matter" will be skipped
%
%   Returns:
%       xyzsInf =           nx3 double, XYZ positions of gray matter electrodes in inflated gifti space. Electrodes in
%                               white matter, contralateral hemisphere, or beyond <distThresh> are returned as [NaN, NaN, NaN] rows.
%
%   HH 2022
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Calculates inflated positions from noninflated pial surface
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
function xyzsInf = getXyzsInfSimple(xyzs, hemi, gii, giiInf, distThresh, destrieuxLabels)

    if nargin < 5 || isempty(distThresh), distThresh = 3; end % maximum distance allowed between valid electrode and gifti vertex
    
    hemi = lower(hemi);
    assert(ismember(hemi, {'r', 'l'}), "hemi must be given as 'r' or 'l'");
    assert(size(gii.vertices, 1) == size(giiInf.vertices, 1), 'gii and giiInf do not have the same number of vertices');

    xyzsInf = nan(size(xyzs));
    for ii = 1:size(xyzs, 1)
        
        if exist('destrieuxLabels', 'var') && ~(startsWith(destrieuxLabels{ii}, sprintf('%sh', hemi)) || contains(destrieuxLabels{ii}, 'Cerebral_White_Matter')), continue; end % must start with lh or rh to be a cortical site
        
        dists = vecnorm(gii.vertices - xyzs(ii, :), 2, 2); % distance from each vertex in gifti to current electrode
        [minDist, idx] = min(dists);
        
        if minDist <= distThresh
            xyzsInf(ii, :) = giiInf.vertices(idx, :);
        end
        
    end
    
end
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: helper function
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
% Edited by HH from color_add_popout_prescaled
%
% Doc for original code:
% function color_add_popout(locs,elcol,msize,th,phi)
% pops electrodes out in th, phi polar coords
% assumes 0<=wts<=1
% all wts should be >=0, if <0, it sets to 0
% kjm 06/20
%
function color_add_custom(locs, wts, cmax, cscale, clip, szrng, marker)

if nargin < 7
    marker = 's';
end

if nargin < 6
    szrng = [2, 25];
end

if nargin < 5 % clip values above clip
    clip = max(wts);
end

if nargin < 4
    cscale = 0.8; % sets the minimum color relative to white (1 = white, 0 = maxcol)
end
assert(numel(szrng)==2 && szrng(2)>szrng(1), 'szrng must be length 2 [minsize, maxsize]');

if length(wts) == 1, wts = repmat(wts, [size(locs, 1), 1]); end % to handle single weight input

wts = wts/clip;
wts(wts<0) = 0; % no negatives
wts(wts>1) = 1; % remove clipped values

if size(cmax, 1)==1
    cmax = repmat(cmax, [size(locs,1), 1]);
end

hold on
for k=1:size(locs,1)
    
    if isnan(wts(k)) || any(isnan(locs(k, :))), continue; end % weight or location missing
    
    mincol = cmax(k, :) + ([1 1 1] - cmax(k, :))*cscale; % min color is between max color and white
    k_col = mincol - (mincol - cmax(k, :))*(wts(k));
    
    if cmax(k, 1) == cmax(k, 2) && cmax(k, 2) == cmax(k, 3) && cmax(k, 1) > 0.15
        edgecol = [0.99, 0.99, 0.99]; % white
        linw = 1.2;
    else
        edgecol = cmax(k, :);
        linw = 0.8;
    end
    
    plot3(locs(k,1), locs(k,2), locs(k,3), marker,...
    'MarkerSize',(szrng(2)-szrng(1))*wts(k)+szrng(1),...
    'LineWidth', linw, ...
    'MarkerEdgeColor', edgecol, ...
    'MarkerFaceColor', k_col)
end


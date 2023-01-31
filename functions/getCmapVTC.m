%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Generates colormap for BPCs, BPSs, and anatomical regions
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
function cm = getCmapVTC(set)
    
    % For BPC curves
    if strcmpi(set, 'bpc')
        
        cm = [0.1, 0,     0.969; % blue
              0.129, 0.62,  0.737; % blue green
              0.466, 0.674, 0.188; % matlab default-5 (green)
              0.694, 0,     0.910]; % purple
         % does not include the gray for unassigned pairs
         
    elseif strcmpi(set, 'bps') % for BPS
        cm = [0.749, 0.345, 0.326;
              0.97, 0.73, 0.39];
              
    
    % For hippocampus, amygdala, VTC, lateral temporal lobe, insula, 
    elseif strcmpi(set, 'anat')
        cm = [0.616, 0.008, 0.031; % dark red, hippocampus
              0.988, 0.749, 0.286; % maximum yellow red, amygdala
              0.992, 1,     0.714; % lemon yellow crayola, ventral temporal
              1,     0.804, 0.698; % apricot, lateral temporal
              0.5,   0.333, 0.222; % coffee, insula
              0.8,   0.8,   0.8 ]; % gray, other
    end
end
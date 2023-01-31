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
function cm = getCmapSpec()

    % Make a nice colormap
    cm1 = repmat([0 0 0], 100, 1);
    cm1(1:40, 1) = 0.7;
    cm1(1:40, 2) = (0.7:-0.6/39:0.1)';
    cm1(1:40, 3) = (0.7:-0.7/39:0)';
    
    cm1(40:100,1) = (0.7:(1-0.7)/60:1)';
    cm1(40:100,2) = (0.1:.9/60:1)';
    
    cm2 = repmat([0 0 0],100,1);
    cm2(1:30,3) = 0.7;
    cm2(1:30,1) = (0.7:-0.7/29:0)';
    cm2(1:30,2) = (0.7:-0.7/29:0)';
    cm2(30:100,3) = (0.7:(1-0.7)/70:1)';
    cm2(30:100,2) = (0:1/70:1)';
    
    cm = [cm2(end:-1:1,:); cm1]; % assemble cm
end
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: wrapper to allow transparent rendering of gifti with SEEG electrodes
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
function renderGiftiElectrodes(electrodesPath, giiPath, viewAngle)
    
    electrodes = readtable(electrodesPath, 'FileType', 'text', 'Delimiter', '\t');
    gii = gifti(giiPath);
    
    figure('Position', [200, 200, 800, 600]);
    ieeg_RenderGifti(gii); alpha 0.3; hold on
    plot3(electrodes.x, electrodes.y, electrodes.z, 'o', 'Color','k', 'MarkerSize',6, 'MarkerFaceColor','w');
    text(electrodes.x, electrodes.y, electrodes.z, electrodes.name, 'Color', 'k');
    hold off
    ieeg_viewLight(viewAngle(1), viewAngle(2));
    
end
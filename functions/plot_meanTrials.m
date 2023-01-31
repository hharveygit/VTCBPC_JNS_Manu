%
%   Plots the mean trial for trials belonging to each stimulation site
%
%   plot_meanTrials(tt, sigdata, sites);
%   plot_meanTrials(tt, sigdata, sites, yspace);
%   plot_meanTrials(tt, sigdata, sites, yspace, varargin);
%
%       tt =            mx1 or 1xm array of timepoints
%       sigdata =       nxm array, each row is a trial and each column corresponds to a timepoint
%       sites =         px2 cell array, where each row corresponds to a site (total p sites). column 1 = name of each site, 
%                       column 2 = row indices in sigdata that belong to that site.
%                       E.g. {'LA1-LA2'} {[1;2;3;4;5;6;7;8;9;10;11;12]}
%                            {'LA2-LA3'} {[13;14;15;16;17;18]}
%                       or px1 cell array, where it assumed that the second column matches the indices of the first
%                       column (no mean)
%                            ...
%       yspace =        (optional) vertical distance between mean trials in the plot. Also length of the scale. 500 by default
%       varargin =      variable input args to plot
%
%   Returns:
%       meanTrials =    pxm array, where each row is a separate mean trial of m timepoints
%                           meanTrials is plotted to new figure
%
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
function [f, meanTrials] = plot_meanTrials(tt, sigdata, sites, yspace, varargin)
    
    if nargin < 4 || isempty(yspace), yspace = 500; end
    tt = tt(:)';

    meanTrials = nan([size(sites, 1), size(sigdata, 2)]);
    
    f = figure('Position', [800, 200, 600, size(sites, 1)*25 + 150]); hold on
    for i = 1:size(meanTrials, 1)
        if size(sites, 2) == 1
            assert(size(sites, 1) == size(sigdata, 1), 'If no indices are given, sigdata must match sites in length');
            meanTrials(i, :) = sigdata(i, :);
        else
            meanTrials(i, :) = mean(sigdata(sites{i, 2}, :), 1);
        end
        adjusted = meanTrials(i, :) - (i-1)*yspace;
        plot(tt, adjusted, varargin{:});
        yline(-(i-1)*yspace, 'Color', 0.5*[1 1 1]);
    end
    
    ylim([-(i+3)*yspace, yspace*4]);
    xline(0, 'Color', 'r');
    set(gca, 'YTick', -yspace*(size(sites, 1)-1:-1:0), 'YTickLabel', flip(sites(:, 1)), 'FontSize', 10);
    xlabel('time (s)');
    
    %plot([0.05 0.05], [yspace*0.5 yspace*1.5], 'k', 'LineWidth', 2);
    %vtext(0.055, yspace, sprintf('%d %sV', yspace, native2unicode(181,'latin1')));
end
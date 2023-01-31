%% Weights input signals by exponential function of time constant tau
%   V = weightExponential(V, tt, tau, invert)
%   V = weightExponential(V, tt, tau)
%       V =         txn, rows are time points, columns are signals
%       tt =        tx1, time points corresponding to rows in V
%       tau =       1x1, time constant for exponential
%       invert =    (optional), true for positive exponential, false for decay (false by default)
%
% HH 2021
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
function V = weightExponential(V, tt, tau, invert)
    if nargin < 4, invert = false; end
    tt = tt(:); % ensure column vector
    assert(size(V, 1) == length(tt), 'V and tt must have the same number of samples');
    
    if invert
        weights = exp(tt./tau);
    else
        weights = exp(-tt./tau);
    end
    
    V = V.*weights;
    
    if invert
        V = V./vecnorm(V); % normalize vectors
        fprintf('normalizing unweighted vectors\n');
    end
end
%
%   V: samples by time
%   rng: [start stop] or [start:stop] indices to remove artifact on
%   method: 'csaps', 'spline' or 'crowther'
%
%   Returns
%       V_clean: same size as V
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
function V_clean = remove_artifact(V, rng, method)
    if nargin < 3
        method = 'csaps';
    end
    
    if numel(rng) == 2 % allow user to input [start end] or [start:end]
        rng = rng(1):rng(end);
    elseif ~(numel(rng) > 2 && isequal(rng, rng(1):rng(end))) % rng is unchanged if rng is start:end vector
        error('rng must be [start stop] or [start:stop] indices of the artifact')
    end
    
    if size(V, 2) > size(V, 1)
        warning('Input V has rows longer than columns -- rows should be time points');
    end
    
    switch method
        case 'csaps'
            V_clean = remove_artifact_csaps(V, rng, 1e-2);
        case 'spline' % non-smoothed cubic spline, deprecated just use csaps with p=1
            V_clean = remove_artifact_spline(V, rng);
        case 'crowther'
            V_clean = remove_artifact_crowther(V, rng);
        otherwise
            error('Method (3rd arg) must be ''csaps'', ''spline'' or ''crowther''');
    end
end

function V_clean = remove_artifact_csaps(V, ix, p, samps) % Replaces artifact with cubic spline interpolation of <samps> # samples before/after
    
    if nargin < 4
        samps = [10, 10];
    end
    
    V_clean = V;
    
    temp_ix = (ix(1)-samps(1)):(ix(end)+samps(2)); % interpolate on 10 samples (~5ms) before & after artifact
    for i=1:size(V_clean, 2)
        temp = V_clean(temp_ix, i);
        
        if all(isnan(temp)), continue; end % ignore nan-rows
        
        ixrng = 1:length(temp_ix); % temp_ix renumerated to begin at 1 (to standardize behavior of smoothing p)
        temp(samps(1)+1:end-samps(2)) = []; % remove artifact
        ixrng(samps(1)+1:end-samps(2)) = [];
        
        vals = csaps(ixrng, temp, p, floor(samps(1)/2+1):ceil(length(temp_ix)-samps(2)/2)); % interpolated values on interval of artifact
        %vals = csaps(ixrng, temp, p, 1:length(temp_ix));
        
        %V_clean(i, ix) = vals; % artifact range only
        V_clean(temp_ix(floor(samps(1)/2+1:ceil(length(temp_ix)-samps(2)/2))), i) = vals; % end halfway through samps
    end

end

function V_clean = remove_artifact_spline(V, ix, samps) % Replaces artifact with cubic spline interpolation of <samps> # samples before/after
    
    warning('Recommend using ''csaps'' instead of ''spline''');
    
    if nargin<3
        samps = 10;
    end
    
    V_clean = V;
    V_clean(ix, :) = nan;
    
    temp_ix = (ix(1)-samps):(ix(end)+samps); % interpolate on 10 samples (~5ms) before & after artifact
    for i=1:size(V_clean, 2)
        temp = V_clean(temp_ix, i);
        y = fillmissing(temp, 'spline', 'SamplePoints', temp_ix);
        V_clean(temp_ix, i) = y;
    end

end

function V_clean = remove_artifact_crowther(V, ix) % ix must be start:stop vector of indices
    
    npoints = length(ix);

    ix_back = ix(end:-1:1) - npoints; % reversed
    ix_forward = ix(end:-1:1) + npoints;
    
    ratios = linspace(1/(npoints+1), (1 - 1/(npoints+1)), npoints)'; % values at k + npoints-k add to 1
    
    V_clean = V;
    V_clean(ix, :) = V(ix_back, :).*ratios(end:-1:1) + V(ix_forward, :).*ratios;
end
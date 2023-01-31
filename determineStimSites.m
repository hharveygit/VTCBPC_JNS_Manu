%% Calculates number of stimulation sites per subject
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
%% Configure subjects

clear;

% Channel names
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6', 'RC6', 'RC7'}, {'RA9', 'RA10'}};
hemis = 'lrlrr'; % hemispheres corresponding to subjects 1-5
subs = struct();
for ii = 1:5
    subs(ii).name = num2str(ii); % de-ID sub names
    subs(ii).chs = subChs{ii};
    subs(ii).hemi = hemis(ii);
end

%% Iterate across subjects to determine number of stim sites

nsites = nan(5, 1);
for subInd = 1:5
    
    %%
    sub = subs(subInd).name;
    chs = subs(subInd).chs;

    subDir = fullfile('data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg');

    mefPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_ieeg.mefd', sub));
    channelsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_channels.tsv', sub));
    eventsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_events.tsv', sub));

    mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
    mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'});

    events = mefObj.evts;

    % Remove all events that don't fit 'eln-el(n+1)' naming
    pairNum = zeros(size(events, 1), 2);
    for kk = 1:length(pairNum)
        pairNum(kk, :) = str2double(regexp(events.electrical_stimulation_site{kk}, '\d*', 'match'));
    end
    events(diff(pairNum, 1, 2) ~= 1, :) = [];

    % Remove all events that correspond to seizure onset zones
    electrodes = readtableRmHyphens(fullfile(subDir, sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub)));
    elecsSoz = electrodes.name(contains(electrodes.seizure_zone, 'SOZ'));
    eventsSoz = any(ismember(split(events.electrical_stimulation_site, '-'), elecsSoz), 2);
    fprintf('Removed %d events at %d sites in SOZ\n', sum(eventsSoz), length(unique(events.electrical_stimulation_site(eventsSoz))));
    events(eventsSoz, :) = [];
    
    sites = groupby(events.electrical_stimulation_site);
    
    % resolve only 4mA or 6mA trials for stim sites with both
    for ii = 1:size(sites, 1)
        idxes = sites{ii, 2};
        trialsCurr = events(idxes, :);
        if length(unique(trialsCurr.electrical_stimulation_current)) == 1, continue; end
        
        trials6ma = idxes(strcmp(trialsCurr.electrical_stimulation_current, '6.0 mA')); % trials at 6 mA stim
        if length(trials6ma) >= 8
            sites{ii, 2} = trials6ma;
            continue
        else % assign the 4.0 mA trials
            sites{ii, 2} = idxes(strcmp(trialsCurr.electrical_stimulation_current, '4.0 mA'));
        end
    end
    sites(cellfun(@length, sites(:, 2)) < 8, :) = []; % exclude sites with < 8 trials
    
    nsites(subInd) = length(sites); % all good stim sites per subject, including the recording electrodes
    
end

fprintf('\n');
for subInd = 1:5
    fprintf('Sub-%d: %d stim sites\n', subInd, nsites(subInd));
end

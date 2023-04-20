%% Plots inflated pial rendering and subject T1 MRI with subject BPC-labeled stim sites for sub1 to include in Figure 2
%
%    If this code is used in a publication, please cite the manuscript:
%    "Electrical stimulation of temporal, limbic circuitry produces multiple
%    distinct responses in human ventral temporal cortex"
%    by H Huang, NM Gregg, G Ojeda Valencia, BH Brinkmann, BN Lundstrom,
%    GA Worrell, KJ Miller, and D Hermes.
%
%    VTCBPC manuscript package: Plots stim sites to subject inflated pial surfaces, colored by consensus BPC
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
%% Configure user parameters, paths, Load data

cm = getCmapVTC('bpc');
cm(5, :) = [0.3, 0.3, 0.3]; % add color for excluded consensus BPCs

sub = '1';
ch = 'LC6';
hemi = 'l';

elecs = readtableRmHyphens(fullfile('data', 'sub-1', 'ses-ieeg01', 'ieeg', 'sub-1_ses-ieeg01_electrodes.tsv'));
gii = gifti(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', sub), sprintf('pial.%s.surf.gii', upper(hemi))));
giiInf = gifti(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', sub), sprintf('inflated.%s.surf.gii', upper(hemi))));
sulc = read_curv(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', sub), 'surf', sprintf('%sh.sulc', lower(hemi))));

xyzsInf = getXyzsInf(elecs, hemi, gii, giiInf, 6);
th = 30*(hemi-111); % angle to view gifti

snrTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', sub), sprintf('sub-%s_ch-%s_BPCSNR.tsv', sub, ch)), 'electrical_stimulation_site', 1);

    
%% Plot and save inflated gifti

recElec = xyzsInf(strcmp(elecs.name, ch), :); % index of recording electrode
        
pairXyzs = ieeg_getPairXyzs(split(snrTable.electrical_stimulation_site, '-', 2), elecs);
pairXyzsInf = getXyzsInfSimple(pairXyzs, hemi, gii, giiInf, 6, snrTable.labelDestrieux); % does not plot white matter electrodes or electrodes in wrong hemisphere

views = [th, 0; -th, 0; th, -90]; % plot 3 views (lateral; medial, inferior); lateral and inferior shown in figure 8
for vv = 1:3

    f = figure('Position', [1000, 100, 500, 350]); % normal gifti
    ieeg_RenderGifti(giiInf, sulc); hold on

    ieeg_viewLight(views(vv, 1), views(vv, 2));
    if vv == 3
        recElecPop = els_popout(recElec, views(vv, 1), views(vv, 2), 30); % pop out a lot when viewing inferior
    else
        recElecPop = els_popout(recElec, views(vv, 1), views(vv, 2), 1); % pop out a little
    end
    plot3(recElecPop(1), recElecPop(2), recElecPop(3), 'o', 'Color', 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

    % color on giftis with consensus BP
    pairXyzsInfPop = els_popout(pairXyzsInf, views(vv, 1), views(vv, 2));
    for bb = 1:max(snrTable.BPC)
        color_add_custom(pairXyzsInfPop(snrTable.BPC == bb, :), ...
                         snrTable.SNR(snrTable.BPC == bb), cm(bb, :), 0.6, 3, [6, 20], 's'); % changed cscale to 0.6 from 0.8
    end

    color_add_custom(pairXyzsInfPop(isnan(snrTable.BPC), :), 0.5, 0.1*[1, 1, 1], 0.2, 1, [2, 10], 's'); % changed size, colors
    hold off
    
    saveas(f, fullfile('output', sprintf('sub-%s', sub), 'inflated', sprintf('sub-%s_ch-%s_infl_view%d_subBpcColors.png', sub, ch, vv)));

end

%% (1) Load interactive GUI for specific subject T1 MRI

niiPath = fullfile('data', 'derivatives', 'freesurfer', 'sub-1', 'sub-1_ses-mri01_T1w_acpc_deFaced.nii');

opts.plotnames = false;
opts.elecWidth = 8;
opts.initialPos = 0;
opts.clim = [0, 0.6];

hands = sliceGUI(niiPath, elecs, opts);

return % pause here, then
% Manually click to Z=-15.0, X=-23.0 on the coronal slice

%% (2) Add BPC stuff for the desired recording channel

% non-assigned stim sites
optsNan.marker = 's';
optsNan.szrng = [3, 7];
addSliceStims(hands, pairXyzs(isnan(snrTable.SNR), :), [], [0.1, 0.1, 0.1], optsNan)

optsStims.clip = 3;
optsStims.marker = 's';
optsStims.szrng = [6, 12];
optsStims.cscale = 0.6;
addSliceStims(hands, pairXyzs, [snrTable.SNR, snrTable.BPC], cm, optsStims);

% Save
kjm_printfig(fullfile('output', sprintf('sub-%s', sub), 'hippAmygSlices', sprintf('sub-%s_ch-%s_slice_subBpcColors', sub, ch)), [30, 10]);

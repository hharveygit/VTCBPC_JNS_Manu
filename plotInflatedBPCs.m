%% Plots stim sites to subject inflated pial surfaces, colored by consensus BPC
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
%% Configure subjects and colormap

clear;
subNames = {'1', '2', '3', '4'}; % Collateral sulcus subjects
subChs = {{'LC6', 'LB6', 'LB7', 'LC5'}, {'RC6'}, {'LB6', 'LB7', 'LC4', 'LC5', 'LD4'}, {'RB5', 'RB6'}};
hemis = 'lrlr';

cm = getCmapVTC('bpc');
cm(5, :) = [0.3, 0.3, 0.3]; % add color for excluded consensus BPCs

%% Iterate over subjects (go through manually if wish to permute BPC order)

for ii = 1:length(subNames)
    elecs = readtableRmHyphens(fullfile('data', sprintf('sub-%s', subNames{ii}), 'ses-ieeg01', 'ieeg', sprintf('sub-%s_ses-ieeg01_electrodes.tsv', subNames{ii})));
    
    gii = gifti(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', subNames{ii}), sprintf('pial.%s.surf.gii', upper(hemis(ii)))));
    giiInf = gifti(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', subNames{ii}), sprintf('inflated.%s.surf.gii', upper(hemis(ii)))));
    sulc = read_curv(fullfile('data', 'derivatives', 'freesurfer', sprintf('sub-%s', subNames{ii}), 'surf', sprintf('%sh.sulc', lower(hemis(ii)))));
    
    xyzsInf = getXyzsInf(elecs, hemis(ii), gii, giiInf, 6);
    th = 30*(hemis(ii)-111); % angle to view gifti
    
    mkdir(fullfile('output', sprintf('sub-%s', subNames{ii}), 'inflated'));
    
    %% Iterate over electrodes
    
    for jj = 1:length(subChs{ii})
        %% Plot and save BPCs
        
        % load SNR values and BPCs
        snrTable = readtableRmHyphens(fullfile('output', sprintf('sub-%s', subNames{ii}), sprintf('sub-%s_ch-%s_BPCSNR.tsv', subNames{ii}, subChs{ii}{jj})), 'electrical_stimulation_site', 1);
        load(fullfile('output', sprintf('sub-%s', subNames{ii}), sprintf('sub-%s_ch-%s_BPCs.mat', subNames{ii}, subChs{ii}{jj})))
        
        % matching subject level to consensus BPCs
        sub2cons = unique([snrTable.BPC(~isnan(snrTable.BPC)), snrTable.consBPC(~isnan(snrTable.BPC))], 'rows');
        sub2cons = sortrows(sub2cons, 2); % sort by consensus BPC to determine shuffling order
        order = sub2cons(:, 1);
        cmConsensus = cm(sub2cons(:, 2), :); % reconstructed colormap to use to account for duplicates
        
        % For color purposes, pad to 4 BPCs
        curvesFilled = zeros(size(curves, 1), 4);
        curvesFilled(:, 1:size(curves, 2)) = curves;
        
        % order by consensus assignments
        bpcsOrdered = curvesFilled(:, order);
        
        % save plots on inflated brain
        f = figure('Position',[100 100 800 600]); hold on
        xlabel('time (s)');
        ylabel('V (unit-normalized)');
        set(gca, 'YTicklabel', []);
        plotTrials(ttBPC, weightExponential(bpcsOrdered, ttBPC, 0.1, true), 0.18, [], cmConsensus, 'LineWidth', 1.5);
        saveas(f, fullfile('output', sprintf('sub-%s', subNames{ii}), 'inflated', sprintf('sub-%s_ch-%s_BPCs_infl.svg', subNames{ii}, subChs{ii}{jj})));
        close(f);
        
        %% Plot and save giftis
        
        % index of recording electrode
        recElec = xyzsInf(strcmp(elecs.name, subChs{ii}{jj}), :); 
        
        % Calculate stim site coordinates on inflated brain
        pairXyzs = ieeg_getPairXyzs(split(snrTable.electrical_stimulation_site, '-', 2), elecs);
        pairXyzsInf = getXyzsInfSimple(pairXyzs, hemis(ii), gii, giiInf, 6, snrTable.labelDestrieux); % does not plot white matter electrodes or electrodes in wrong hemisphere
        
        views = [th, 0; -th, 0; th, -90]; % plot 3 views (lateral; medial, inferior); lateral and inferior views shown in figure 8
        for vv = 1:3
            
            f = figure('Position', [1000, 100, 400, 300]); % normal gifti
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
            for bb = 1:5
                color_add_custom(pairXyzsInfPop(snrTable.consBPC == bb, :), ...
                                 snrTable.SNR(snrTable.consBPC == bb), cm(bb, :), 0.6, 3, [6, 20], 's'); % changed cscale to 0.6 from 0.8
            end
            
            color_add_custom(pairXyzsInfPop(isnan(snrTable.consBPC), :), 0.5, 0.1*[1, 1, 1], 0.2, 1, [2, 10], 's'); % changed size, colors
            hold off

            saveas(f, fullfile('output', sprintf('sub-%s', subNames{ii}), 'inflated', sprintf('sub-%s_ch-%s_infl_view%d.png', subNames{ii}, subChs{ii}{jj}, vv)));
            close(f);
            fprintf('Stim sites plotted and saved to inflated brain for subject %s, ch %s\n', subNames{ii}, subChs{ii}{jj});
            
        end
        
    end
    
end
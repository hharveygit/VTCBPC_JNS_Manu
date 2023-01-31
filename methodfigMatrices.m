%% To use for saving projection, significance matrices and NNMF outputs for Figure 1. Uses data from sub1, electrode 1
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
%% Configure paths, load mef data

clear;
set(0, 'DefaultFigureRenderer', 'painters');

mkdir('output/methodsFigure');

sub = '1';
ch = 'LC6';

% location of subject data. Mef needs absolute path, hence pwd
subDir = fullfile(pwd, 'data', sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg');

mefPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_ieeg.mefd', sub));
channelsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_channels.tsv', sub));
eventsPath = fullfile(subDir, sprintf('sub-%s_ses-ieeg01_task-ccep_run-01_events.tsv', sub));

mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'});
mefObj.loadMefTrials([-1, 2]);
mefObj.car(true); % by 64-block
mefObj.pruneChannels(ch);
mefObj.removeLN('SpectrumEstimation');
mefObj.subtractBaseline([-0.5, -0.05], 'median');
mefObj.plotInputs([], [], 200);

events = mefObj.evts;
tt = mefObj.tt;
srate = mefObj.srate;
data = mefObj.data;

% Remove all events that don't fit 'eln-el(n+1)'
pairNum = zeros(size(events, 1), 2);
for kk = 1:length(pairNum)
    pairNum(kk, :) = str2double(regexp(events.electrical_stimulation_site{kk}, '\d*', 'match'));
end
data(:, :, diff(pairNum, 1, 2) ~= 1) = [];
events(diff(pairNum, 1, 2) ~= 1, :) = [];

% Remove all events that correspond to seizure onset zones
electrodes = readtableRmHyphens(fullfile(subDir, sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub)));
elecsSoz = electrodes.name(contains(electrodes.seizure_zone, 'SOZ'));
eventsSoz = any(ismember(split(events.electrical_stimulation_site, '-'), elecsSoz), 2);
fprintf('Removed %d events at %d sites in SOZ\n', sum(eventsSoz), length(unique(events.electrical_stimulation_site(eventsSoz))));
data(:, :, eventsSoz) = [];
events(eventsSoz, :) = [];

assert(size(data, 3) == height(events), 'data - events mismatch');

% Pull out data from channel
fprintf('Getting data from channel %s\n', ch);

sigdata = squeeze(data(strcmp(mefObj.channels.name, ch), :, :))';
events_ch = events;

% remove stimulated channels from current events
sigdata(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
events_ch(any(strcmp(ch, split(events.electrical_stimulation_site, '-')), 2), :) = [];
sites = groupby(events_ch.electrical_stimulation_site);

% resolve only 4mA or 6mA trials for stim sites with both
for ii = 1:size(sites, 1)
    idxes = sites{ii, 2};
    trialsCurr = events_ch(idxes, :);
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

%% Create projection matrix

disp('Getting BPCs');
tau = 0.1;

pairTypes = struct('pair', sites(:, 1), 'indices', sites(:, 2));
V = sigdata(:, tt >= 0.011 & tt < 0.5)';
ttBPC = tt(tt >= 0.011 & tt < 0.5);
V = weightExponential(V, tt(tt >= 0.011 & tt < 0.5), tau); % weight by decaying exponential of time constant 100ms

Vnorm = V./vecnorm(V, 2, 1);
P = V'*Vnorm;
P(1:length(P)+1:end) = nan;
figure;
imagescNaN(1:length(P), 1:length(P), P);
colormap(parula);
caxis([-2000, 2000]);
axis square

% Add crosshair to mark trials belonging to a stim site pair
yline(pairTypes(14).indices(1), 'Color', 'r'); yline(pairTypes(14).indices(end), 'Color', 'r');
xline(pairTypes(18).indices(1), 'Color', 'r'); xline(pairTypes(18).indices(end), 'Color', 'r');

saveas(gcf, fullfile('output', 'methodsFigure', 'projectionmatrix'), 'svg');
saveas(gcf, fullfile('output', 'methodsFigure', 'projectionmatrix'), 'png');

%% Create significance matrix

Xi = nan(length(pairTypes));
for ii = 1:length(pairTypes)
    for jj = 1:length(pairTypes)
        projs = P(pairTypes(ii).indices, pairTypes(jj).indices);
        projs = projs(:); % flatten
        projs(isnan(projs)) = [];
        
        Xi(ii, jj) = mean(projs) / (std(projs)/sqrt(length(projs))); % calculate t value
    end
end

Xi(Xi < 0) = 0;

figure;
imagescNaN(1:length(Xi), 1:length(Xi), Xi);
colormap(flipud(gray));
caxis([0, 30]);
axis square

saveas(gcf, fullfile('output', 'methodsFigure', 'significancematrix'), 'svg');
saveas(gcf, fullfile('output', 'methodsFigure', 'significancematrix'), 'png');

%% Create NNMF decomposition

% multiple run-throughs of  nnmf
tmp_mat = []; tmp_err = [];
for ii = 1:100
    [tmp_mat(ii).W, tmp_mat(ii).H, tmp_err(ii)] = kjm_nnmf(Xi, 3);
end

[~, ind] = min(tmp_err); % find min error

W = tmp_mat(ind).W;
H = tmp_mat(ind).H;

figure;
imagesc(W); colormap(flipud(gray));
caxis([0, prctile(W(:), 95)]);
axis equal; axis off
saveas(gcf, fullfile('output', 'methodsFigure', 'W'), 'svg');
saveas(gcf, fullfile('output', 'methodsFigure', 'W'), 'png');

figure;
imagesc(H); colormap(flipud(gray));
caxis([0, prctile(H(:), 95)]);
axis equal; axis off
saveas(gcf, fullfile('output', 'methodsFigure', 'H'), 'svg');
saveas(gcf, fullfile('output', 'methodsFigure', 'H'), 'png');

Hthresh = zeros(size(H));
thresh = 1/(2*sqrt(length(pairTypes)));
for ii = 1:size(H, 2)
    [Hmax, ind] = max(H(:, ii));
    Hthresh(ind, ii) = Hmax;
end
Hthresh(Hthresh < thresh) = 0; % do not assign any BPC if less than threshold
Hthresh(Hthresh > 0) = 1; % binarize the rest

figure;
imagesc(Hthresh); colormap(flipud(gray));
%caxis([0, prctile(H(:), 95)]);
axis equal; axis off
saveas(gcf, fullfile('output', 'methodsFigure', 'Hthresh'), 'svg');
saveas(gcf, fullfile('output', 'methodsFigure', 'Hthresh'), 'png');



%%
function [W,H, newerr] = kjm_nnmf(V,rdim)
% kjm_nnmf - non-negative matrix factorization 
% Copied from bpc_identify as subfunction
% 
% SYNTAX:
% [W,H] = kjm_nmfsc(V, rdim);
%
% INPUTS:
% V          - data matrix (N x T)
% rdim       - (M) number of 'reduced dimension' components (inner dimension of factorization)
% 
% OUTPUTS:
% W           - first non-negative factor in V~WH decomposition (N x rdim)
% H           - second non-negative factor in V~WH decomposition (rdim x T)
% newerr      - latest estimate of error of in estimation
%
% file is based upon elements of file nmfsc by Patrik Hoyer at https://github.com/aludnam/MATLAB/tree/master/nmfpack
% ref is Hoyer, 2004, 'Non-negative matrix factorization with sparseness constraints' Journal of Machine Learning Research  5:1457-1469, 2004.
% kjm 5/2020


%% defaults and data check

    conv_thresh=1e-5; %convergence threshold

    % Check that we have non-negative data
    if min(V(:))<0, error('Negative values in data!'); end

    % Globally rescale data to avoid potential overflow/underflow
    V = V/max(V(:));
    
%% Create initial random positive matrices
    W = abs(randn(size(V,1),rdim)); 
    H = abs(randn(rdim,size(V,2)));
    
    % normalizes so each row of H has unit L2 norm
    H = H./(sqrt(sum(H.^2,2))*ones(1,size(V,2)));

%% initialize iterative process
    % Calculate initial error (E)
    errhistory = 0.5*sum(sum((V-W*H).^2));
    loop_continue=1;
    loop_num=0;
    
%% iterative convergence
while loop_continue==1

    loop_num=loop_num+1;
        
    %% ----- Update H ---------------------------------------
        % Update using standard NMF multiplicative update rule, (element-wise multiplication of ratios of positive numbers ensures positivity)
        H = H.*(W.'*V)./(W.'*W*H + 1e-9); %adds 1e-9 so doesn't explode if divide by zero

        % Renormalize so rows of H have unit norm
        Hnorms = sqrt(sum(H'.^2));
        H = H./(Hnorms'*ones(1,size(V,2)));
        W = W.*(ones(size(V,1),1)*Hnorms); % pushes scaling from rows of H onto columns of W
	    
    %% ----- Update W ---------------------------------------
        % Update using standard NNMF multiplicative update rule 
        % (element-wise multiplication of ratios of positive numbers ensures positivity)
        W = W.*(V*H.')./(W*H*H.' + 1e-9);	
	        
    %% ----- Calculate error (E) ---------------------------------------
        newerr = norm(V-W*H, 'fro')^2;
        errhistory = [errhistory newerr];

        steperrratio=abs(diff(errhistory([-1 0]+end)))/newerr;
        if steperrratio<conv_thresh, 
%             disp(['ratio of change in error to error below threshold of ' num2str(conv_thresh)]);
            loop_continue=0;
        end
%     
%         if or(mod(loop_num,1000)==1,loop_continue==0)
%             disp(['Iteration: ' num2str(loop_num)...
%                 ', error term = ' num2str(newerr)...
%                 ', error step ratio = ' num2str(steperrratio)])
%         end
end

[tmp,pp]=sort(sum(W.^2,1),'descend');
W=W(:,pp);H=H(pp,:);
end
   
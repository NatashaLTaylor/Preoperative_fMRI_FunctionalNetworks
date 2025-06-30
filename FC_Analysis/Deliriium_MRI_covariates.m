%% check the MRI scanner covariate in delirium fMRI dataset

%load in
load('updated_fc_results.mat');
%import "mri_scanner_current.csv" file - 

%avg FC across each region 
fc_edges_all;

%take lower tril of the subjects fc matrix
fc_edges_all_tril = zeros(120,502,502);
for ii=1:120
    fc_edges_all_tril(ii,:,:) = tril(squeeze(fc_all(ii,:,:)));
end
%replace 0 with NaN
fc_edges_all_tril(fc_edges_all_tril==0)=NaN;
%calculate new avg. ignoring NaN (these are repeated fc)
fc_avg_roi2 = squeeze(nanmean(fc_edges_all_tril,2));


fc_avg_roi = mean(fc_all,2);
fc_avg_roi = squeeze(fc_avg_roi); %120 x 502

%sig difference 
p = zeros(502, 1);
%stats = struct[502,:]; % Initialize an empty structure

for i = 1:502
    [p(i), ~, stats(i)] = anova1(fc_avg_roi(:, i), mriscannercurrent.mri_scanner_used, 'off'); % Suppress ANOVA plot
end

for i = 1:502
    [p(i), ~, stats(i)] = anova1(fc_avg_roi2(:, i), mriscannercurrent.mri_scanner_used, 'off'); % Suppress ANOVA plot
end
%determine significantly different ROIs

loc_sig = (find(p<0.05));
sig_loc_rois = zeros(502,1);
sig_loc_rois(loc_sig,1)=1;

for i=1:502
    [p_kw(i,:), ~, stats_kw(i)] = kruskalwallis(fc_avg_roi(:, i), mriscannercurrent.mri_scanner_used, 'off'); % Suppress plot
end
multcompare(stats_kw(1));

loc_sig_kw = (find(p<0.05));
sig_loc_rois_kw = zeros(502,1);
sig_loc_rois_kw(loc_sig_kw,1)=1;


%% Two-way ANOVA - need to check if the group of scanner significantly determines differences in FC & the delirium outcome

% 120 participants x avg. fc for each ROI
fc_values = fc_avg_roi2;
fc_values = fc_values(:); % Reshape to a column vector (120 * 502 x 1)

%relevant clinical info
mri_scanner = mriscannercurrent.mri_scanner_used;
delirium_outcome = mriscannercurrent.bin_delirium;
%need to be categorical values
mri_scanner = categorical(mri_scanner);
delirium_outcome = categorical(delirium_outcome);
%Repeat Grouping Variables: For 120 participants for each ROI
% Assuming mri_scanner and delirium_outcome are vectors of size 120
mri_scanner = repmat(mri_scanner, 502, 1); % Repeat for each region
delirium_outcome = repmat(delirium_outcome, 502, 1); % Repeat for each region
%Add a Region Factor (Optional): If you want to account for differences between regions in your analysis, you can include a third factor for regions:
regions = repmat((1:502)', 120, 1); % Region index repeated for all participants


%unbalance two-way anova
%Run Two-Way ANOVA: For a two-way ANOVA with MRI Scanner and Delirium Outcome:
[p, tbl, stats] = anovan(fc_values, {mri_scanner, delirium_outcome}, ...
                         'model', 'interaction', ...
                         'varnames', {'MRI Scanner', 'Delirium Outcome'});
multcompare(stats, 'Dimension', 1); % Pairwise comparisons for MRI Scanner

[c,m,h,gnames] = multcompare(stats);
tbl_multcompare = array2table(m,"RowNames",gnames, ...
    "VariableNames",["Mean","Standard Error"])
% check effect size of the groups - partial eta squared
% Extract the variable names for effects (skip header and 'Total' row)
effect_names = tbl(2:end-1, 1); % Exclude first row (header) and last row (Total)

% Convert cell array to proper format for table variable names
effect_names = string(effect_names); % Ensure names are strings or categorical

% Compute partial eta-squared (as before)
ss_effect = cell2mat(tbl(2:end-1, 2)); % Extract SS for effects
ss_error = tbl{end, 2};                % Extract SS for error
eta_squared_partial = ss_effect ./ (ss_effect + ss_error);

% Display the results in a table
disp('Partial Eta-Squared:');
disp(array2table(eta_squared_partial.', 'VariableNames', {'MRI Scanner', 'Delirium Outcome', 'MRI Scanner:Outcome', 'Error'}));


% cohen's d
cohen_f = sqrt(eta_squared_partial ./ (1 - eta_squared_partial));
disp('Cohen''s f:');
disp(array2table(cohen_f.', 'VariableNames',{'MRI Scanner', 'Delirium Outcome', 'MRI Scanner:Outcome', 'Error'}));


%if including regions as a factor
[p_2anova, tbl_2anova, stats_2anova] = anovan(fc_values, {mri_scanner, delirium_outcome, regions}, ...
                         'model', 'interaction', ...
                         'varnames', {'MRI Scanner', 'Delirium Outcome', 'Regions'});
%interaction plto
interactionplot(fc_values, {mri_scanner, delirium_outcome});
multcompare(stats_2anova, 'Dimension', 1); % Pairwise comparisons for MRI Scanner
multcompare(stats_2anova, 'Dimension', 2); % Pairwise comparisons for Delirium Outcome

% Two-way ANOVA - brain regions as a random variable
[p, tbl, stats] = anovan(fc_values, {mri_scanner, delirium_outcome, regions}, ...
                         'model', 'interaction', ...
                         'random', 3, ... % Specify regions as random if necessary
                         'varnames', {'MRI Scanner', 'Delirium Outcome', 'Regions'});

%% Group Diff between MRI scanner 0 and 3

%split them into two groups
grp_3 = mriscannercurrent.mri_scanner_used==3;
grp_mri_3_fc = fc_all(grp_3==1,:,:);

grp_0 = mriscannercurrent.mri_scanner_used==0;
grp_mri_0_fc = fc_all(grp_0==1,:,:);

%determine if significant difference - permutation test
% flatten and loop through
nROIs = size(fc_all,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end

flat_fc_grp_3 = reshape(grp_mri_3_fc,77,nROIs*nROIs);
fc_edges_grp3 = flat_fc_grp_3(:,template);
flat_fc_grp_0 = reshape(grp_mri_0_fc,35,nROIs*nROIs);
fc_edges_grp0 = flat_fc_grp_0(:,template);
nEdges = size(fc_edges_grp3,2);

% run permutation across fc edges - between grp 0 and grp 3 for MRI scanner
% type
for i=1:nEdges
    [sig_fc_edges_mri(i,:),pval_fc_edges_mri(i,:)]=perm_code(fc_edges_grp0(:,i),fc_edges_grp3(:,i),1000);
end


%% run regression from BOLD signal 

%load in time-series that match 
for ii=1:length(filename)
    %get 'sub-...' from the name
    if ii<10
        subject_file= sprintf('%s%d%s','sub-00',ii,'_ses-1_timeseries.mat'); %filename of time-series data
    elseif ii<100
        subject_file= sprintf('%s%d%s','sub-0',ii,'_ses-1_timeseries.mat');
    else
        subject_file= sprintf('%s%d%s','sub-',ii,'_ses-1_timeseries.mat');
    end
    load([subject_file]); %load in time-series; ts variable
    ts = ts(6:end -5,:); %remove first/last 5 time-points for noise
    

    %functional connectivity
    ts_corr = corr(ts); %ROI x ROI matrix
    fc_all(ii,:,:)=ts_corr; %all subjects FC into one
end



subjectV1MRIdata.sub_id()
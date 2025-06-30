%% DMN Analysis


%% 1. FC Analysis 
load('schaef_400_fc_all.mat') %load in 400 network schaef
% load clinical data
load('C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Demographic_Data\V1_MRI_all_demos.mat');
%find IDs that have 'NaN'
sub_remove = find(subjectV1MRIdata.overall_delirious_ever=='NA');
subjectV1MRIdata(sub_remove,:)=[];
bin_delirium_all=table2array(subjectV1MRIdata(1:size(subjectV1MRIdata,1),"bin_delirium"));
delirium_sub = find(bin_delirium_all==1);
health_sub =bin_delirium_all~=1; 
health_sub = find(health_sub==1);
fc_all(sub_remove,:,:)=[];

%remove subjects without peak del. data
fc_all_avg_dmn(sub_remove,:,:)=[];
fc_del = fc_all_avg_dmn(delirium_sub,:,:);
fc_health = fc_all_avg_dmn(health_sub,:,:);
%DMN ids
[149:194,358:390]
fc_del_dmn = fc_del(:,:,[149:194,358:390]);
fc_health_dmn = fc_health(:,:,[149:194,358:390]);
fc_dmn = fc_all(:,:,[149:194,358:390]);

% just look at the difference when ran permutation across entire brain

% fc dmn within
fc_del_within_dmn = fc_del(:,[149:194,358:390],[149:194,358:390]);
sig_fc_del_within_dmn = sig_fc_del([149:194,358:390],[149:194,358:390]);
fc_health_within_dmn = fc_health(:,[149:194,358:390],[149:194,358:390]);
sig_fc_health_within_dmn = sig_fc_health([149:194,358:390],[149:194,358:390]);

diff_fc_del_health_within_dmn = diff_fc_del_health([149:194,358:390],[149:194,358:390]);

%load in DMN order subnets
load('dmn_order_new.mat') %sub-regions of DMN
load('dmn_order_subnets.mat') %the sub-network DMN order


%reorder to 
order_sig_diff_fc_dmn = sig_diff_fc_dmn(order_8_nets(1:400,2),new_order_dmn);
writematrix(order_sig_diff_fc_dmn,'order_diff_dmn_wholebrain.csv')

%sig_diff_fc_dmn
aa = mean(sig_diff_fc_dmn(1:400,:),2);
RB_surf_schaef_pink(aa)


%% Re-run permutation test - 12/09/24
%need to remove self connections
nROIs = size(fc_del_dmn,2);
flat_fc_del_dmn = reshape(fc_del_dmn,size(fc_del_dmn,1),nROIs*79);
flat_fc_health_dmn = reshape(fc_health_dmn,size(fc_health_dmn,1),nROIs*79);

columns_to_remove = all(flat_fc_del_dmn == 1);
flat_fc_del_dmn(:, columns_to_remove) = [];
columns_to_remove = all(flat_fc_health_dmn == 1);
flat_fc_health_dmn(:,columns_to_remove)=[];
locs_self = find(columns_to_remove==1);
%run permutation for flatten fc dmn to whole brain
for i=1:size(flat_fc_health_dmn,2)
    [sig_fc_dmn_brain(i,:),pval_fc_dmn_brain(i,:)]=perm_code(flat_fc_del_dmn(:,i),flat_fc_health_dmn(:,i),1000);
end

sig_fc_del_dmn_flat = sig_fc_dmn_brain.*(squeeze(mean(flat_fc_del_dmn))');
sig_fc_health_dmn_flat = sig_fc_dmn_brain.*(squeeze(mean(flat_fc_health_dmn))');

sig_diff_fc_dmn = sig_fc_del_dmn - sig_fc_health_dmn;

%add back in the removed elements
numNewElements = length(locs_self);
sig_fc_health_dmn_flat2 = zeros(length(sig_fc_health_dmn_flat)+numNewElements,1);
originalIndex =1;
newIndex = 1;
for i = 1:length(sig_fc_health_dmn_flat) + numNewElements
    if ismember(i, locs_self)
        sig_fc_health_dmn_flat2(i) = 0;  % Insert zero at the specified position
    else
        sig_fc_health_dmn_flat2(i) = sig_fc_health_dmn_flat(originalIndex);  % Copy the original element
        originalIndex = originalIndex + 1;
    end
end
sig_fc_health_dmn = reshape(sig_fc_health_dmn_flat2,502,79);

sig_diff_fc_dmn = sig_fc_del_dmn - sig_fc_health_dmn;
%just looking at the difference from the whole brain permutation
sig_diff_fc_dmn2 = diff_fc_del_health(:,[149:194,358:390]);



% take avg. of sig. connections across DMN nodes - do it by group into
%% networks avg. connections to DMN
a = sig_fc_del_dmn;
a(a==0)=NaN; %replaces zeros with Nan
avg_sig_fc_del_dmn = mean(a,2,'omitnan'); %calculate mean, ignoring NaNs
[left, right] = RB_surf_schaef(avg_sig_fc_del_dmn(1:400),'nan_mean_sig_fc_del_dmn')

b = sig_fc_health_dmn;
b(b==0)=NaN; %replaces zeros with Nan
avg_sig_fc_health_dmn = mean(b,2,'omitnan'); %calculate mean, ignoring NaNs
[left, right] = RB_surf_schaef(avg_sig_fc_health_dmn(1:400),'nan_mean_sig_fc_health_dmn')

c = sig_diff_fc_dmn;
c(c==0)=NaN; %replaces zeros with Nan
avg_diff_sig_fc_dmn = mean(c,2,'omitnan'); %calculate mean, ignoring NaNs
[left, right] = RB_surf_schaef(avg_diff_sig_fc_dmn(1:400),'nan_mean_diff_sig_fc_dmn')


%take avg. from the networks


%% DMN avg. connectivity

avg_dmn_del = mean(fc_del_dmn,3); %avg. across all DMN nodes
avg_dmn_health = mean(fc_health_dmn,3); %avg. across all DMN nodes
for i=1:size(avg_dmn_del,2)
    [sig_dmn_avg(i,:),pval_dmn(i,:)]=perm_code(avg_dmn_health(:,i),avg_dmn_health(:,i),1000);
end
sig_avg_dmn_del = sig_dmn_avg'.*mean(avg_dmn_del);
sig_avg_dmn_health = sig_dmn_avg'.*mean(avg_dmn_health); %size 1 x 502

%% Within DMN Analysis

figure
set(gcf,'color','w')
imagesc(sig_diff_fc_dmn([149:194,358:390],:))
title('Sig. Difference DMN Connections Del - Health')
colormap(Color_Salmon_Blue2)
hold on
line([0,79],[45,45],'Color','black','LineWidth',1) %hemisphere line
line([45,45],[79,0],'Color','black','LineWidth',1) %hemisphere line
%additional lines
 hold on
line([0,79],[18,18],'Color','black','LineWidth',1)
line([18,18],[0,79],'Color','black','LineWidth',1)
line([0,79],[39,39],'Color','black','LineWidth',1)
line([39,39],[0,79],'Color','black','LineWidth',1)
line([0,79],[5,5],'Color','black','LineWidth',1)
line([5,5],[0,79],'Color','black','LineWidth',1)
line([0,79],[12,12],'Color','black','LineWidth',1)
line([12,12],[0,79],'Color','black','LineWidth',1)
line([24,24],[0,79],'Color','black','LineWidth',1)
line([0,79],[24,24],'Color','black','LineWidth',1)
line([39,39],[0,79],'Color','black','LineWidth',1)
line([51,51],[0,79],'Color','black','LineWidth',1)
line([0,79],[51,51],'Color','black','LineWidth',1)
line([0,79],[56,56],'Color','black','LineWidth',1)
line([56,56],[0,79],'Color','black','LineWidth',1)
line([0,79],[62,62],'Color','black','LineWidth',1)
line([62,62],[0,79],'Color','black','LineWidth',1)
line([0,79],[66,66],'Color','black','LineWidth',1)
line([66,66],[0,79],'Color','black','LineWidth',1)
line([0,79],[73,73],'Color','black','LineWidth',1)
line([73,73],[0,79],'Color','black','LineWidth',1)


%correlation
fc_del_within_dmn = fc_del_dmn(:,[149:194,358:390],:); %take 79 DMN nodes to itself
fc_health_within_dmn = fc_health_dmn(:,[149:194,358:390],:); %take 79 DMN nodes to itself

fc_dmn_all_within = fc_dmn(:,[149:194,358:390],:);
flat_fc_dmn_all_within = reshape(fc_dmn_all_within,120,79*79);
[pval_permute_within_dmn,orig_corr_within_dmn,orig_pval_within_dmn,null_corr_within_dmn] = perm_1d_corr(peak_drs_overall,flat_fc_dmn_all_within,1000);
sig_pval_permute_corr_within_dmn = double(pval_permute_within_dmn<0.05);
sig_permute_corr_within_dmn = sig_pval_permute_corr_within_dmn.*orig_corr_within_dmn;
reshape_sig_permute_corr_within_dmn = reshape(sig_permute_corr_within_dmn,79,79);

%remove FC edges to self, as this is 1
find_loc = flat_fc_dmn_all_within==1;
a = find(flat_fc_dmn_all_within==1);
flat_flat_fc_dmn_all_within = reshape(flat_fc_dmn_all_within,120*6241,1);
removed_within_fc_dmn = flat_flat_fc_dmn_all_within;
removed_within_fc_dmn(a)=[]; %removed all "to self connections"
reshape_remove_within_fc_dmn = reshape(removed_within_fc_dmn,120,6162);


% only take lower triangle - unique FC edges
nROIs = size(reshape_fc_dmn_within,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
fc_edges_dmn_within = flat_fc_dmn_all_within(:,template); %flat across fc edges only

avg_fc_edges_dmn_within = mean(fc_edges_dmn_within,1);
mat = matfy(avg_fc_edges_dmn_within',79); %

%% load in .csv files from r-studio

permute_coeffs_within_dmn= readmatrix("permute_lm_within_dmn_coeffs_pval.csv"); %first column 1 coeff, 2nd column is 
sig_permute = double(permute_coeffs_within_dmn(:,2)<0.05);
sig_coefs_within_dmn = sig_permute.*permute_coeffs_within_dmn(:,1);
sig_coefs_within_dmn(3082,:)=[]; %as this is the predictions of nsqip-d 
mat_sig_coefs_within_dmn = matify(sig_coefs_within_dmn,79);

figure
set(gcf,'color','w')
imagesc(mat_sig_coefs_within_dmn)
title('Sig. Permute Rho within DMN connections vs Peak DRS')
hold on
line([0,79.5],[45,45],'Color','black','LineWidth',1) %hemisphere line
line([45,45],[79.5,0],'Color','black','LineWidth',1) %hemisphere line
colormap(Color_Salmon_Blue)

figure
set(gcf,'color','w')
imagesc(mat_sig_coefs(new_order_dmn,new_order_dmn))
title('Sig. Permute LM within DMN connections vs Peak DRS')
hold on
%order DMN into regional grps
line([0,502],[9.5,9.5],'Color','black','LineWidth',1) %IPL
line([9.5,9.5],[0,502],'Color','black','LineWidth',1) %IPL
line([0,502],[25.5,25.5],'Color','black','LineWidth',1) %PFCd
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %PFCd
line([0,502],[37.5,37.5],'Color','black','LineWidth',1) %PFCm
line([37.5,37.5],[0,502],'Color','black','LineWidth',1) %PFCm
line([0,502],[47.5,47.5],'Color','black','LineWidth',1) %PFCl/v
line([47.5,47.5],[0,502],'Color','black','LineWidth',1) %PFCl/v
line([0,502],[59.5,59.5],'Color','black','LineWidth',1) %pCUNPCC
line([59.5,59.5],[0,502],'Color','black','LineWidth',1) %pCUNPCC
line([0,502],[69.5,69.5],'Color','black','LineWidth',1) %Temp lobe
line([69.5,69.5],[0,502],'Color','black','LineWidth',1) %Temp lobe
line([0,502],[74.5,74.5],'Color','black','LineWidth',1) %DefaultC_Rsp
line([74.5,74.5],[0,502],'Color','black','LineWidth',1) %DefaultC_Rsp
colormap(Color_Salmon_Blue)
colorbar

figure
set(gcf,'color','w')
imagesc(mat_sig_coefs_within_dmn(order_subnets,order_subnets))
title('Sig. Permute LM within DMN connections vs Peak DRS')
hold on
colormap(Color_Salmon_Blue)
line([0,502],[12.5,12.5],'Color','black','LineWidth',1) %Core subnet
line([12.5,12.5],[0,502],'Color','black','LineWidth',1) %Core subnet
line([0,502],[41.5,41.5],'Color','black','LineWidth',1) %Med Temp subnet
line([41.5,41.5],[0,502],'Color','black','LineWidth',1) %Med Temp subnet
%rest dors/med pFC
colorbar
xticklabels([])
yticklabels([])
% Remove x and y ticks
xticks([]);
yticks([]);

%% Permutation Testing - 5/06/24

%permutation of the correlation - fix the way p-value is calculated in this
for i=1:size(flat_DMN_FC,2)
    [sig_permute_corr(i,:),pval_permute_corr(i,:),orig_delta_permute_corr(i,:)] = permutation_correlate(flat_DMN_FC(:,i),peak_drs_overall,1000);
end


[pval_permute,orig_corr,orig_pval,null_corr] = perm_1d_corr(peak_drs_overall,flat_DMN_FC,1000);
sig_pval_permute_corr = double(pval_permute<0.05);
sig_permute_corr = sig_pval_permute_corr.*orig_corr;
reshape_sig_permute_corr = reshape(sig_permute_corr,502,79);


%% FDR Corrected P-values -
% only needs to be done on the correlated, don't think it's necessary for
% permuted p-values
[FDR,Q,aPrioriProb] = mafdr(PValues,);


%% Extracted files to run permuted linear model
%need to flatten fc_dmn - but remove the correlations to self
a = reshape(fc_dmn,120,502*79);
columns_to_remove = all(a == 1);
a(:, columns_to_remove) = []; %120 x10020
b = 1:size(a,2);
c = vertcat(b,a);

writematrix(c, 'fc_dmn_rest_brain_all_subjects.csv');



%% Permuted Linear Model from R studio, DMN connectivity to rest of brain
permute_coeffs_dmn_fc= readmatrix("permute_lm_fc_dmn_coeffs_pval.csv"); %first column 1 coeff, 2nd column is 
%removed all 79 self connections
sig_permute = double(permute_coeffs_dmn_fc(:,2)<0.05);
sig_coefs = sig_permute.*permute_coeffs_dmn_fc(:,1);
sig_coefs(length(sig_coefs),:)=[]; %as this is the predictions of nsqip-d 
%add in zeros to corrs with self
locs_self = find(columns_to_remove==1);
numNewElements = length(locs_self);
B = zeros(length(sig_coefs)+numNewElements,1);
originalIndex =1;
newIndex = 1;

for i = 1:length(sig_coefs) + numNewElements
    if ismember(i, locs_self)
        B(i) = 0;  % Insert zero at the specified position
    else
        B(i) = sig_coefs(originalIndex);  % Copy the original element
        originalIndex = originalIndex + 1;
    end
end

mat_sig_coefs_lm = reshape(B,502,79); %sig. beta coefs from permute lm dmn fc to brain
sig_coefs_dmn_fc_lm = mat_sig_coefs_lm;
load('schaef_order.mat')
load('dmn_order_new.mat')
load('colormaps2.mat')

figure
set(gcf,'color','w')
imagesc(mat_sig_coefs_lm(voltron_order(:,3),:))
title('Sig. Permuted Coefs FC DMN vs peak drs')
colormap(Color_Salmon_Blue)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)
line([45.5,45.5],[0,502],'Color','black','LineWidth',1) %hemisphere line
%reordered
figure
set(gcf,'color','w')
imagesc(mat_sig_coefs_lm(voltron_order(:,3),new_order_dmn))
title('Sig. Permuted Coefs FC DMN ordered vs peak drs')
colormap(Color_Salmon_Blue)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)
%order DMN into regional grps
%order DMN into regional grps
line([9.5,9.5],[0,502],'Color','black','LineWidth',1) %IPL
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %PFCd
line([37.5,37.5],[0,502],'Color','black','LineWidth',1) %PFCm
line([47.5,47.5],[0,502],'Color','black','LineWidth',1) %PFCl/v
line([59.5,59.5],[0,502],'Color','black','LineWidth',1) %pCUNPCC
line([69.5,69.5],[0,502],'Color','black','LineWidth',1) %Temp lobe
line([74.5,74.5],[0,502],'Color','black','LineWidth',1) %DefaultC_Rsp
% ordered into sub-networks of DMN
figure
set(gcf,'color','w')
imagesc(mat_sig_coefs_lm(voltron_order(:,3),order_subnets))
title('Sig. Permuted Coefs FC DMN ordered vs peak drs')
colormap(Color_Salmon_Blue)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)
line([12.5,12.5],[0,502],'Color','black','LineWidth',1) %Core subnet
line([41.5,41.5],[0,502],'Color','black','LineWidth',1) %Med Temp subnet
%rest dors/med pFC
colorbar
xticklabels([])
yticklabels([])
% Remove x and y ticks
xticks([]);
yticks([]);

figure
set(gcf,'color','w')
xx = mat_sig_coefs_lm';
imagesc(xx(new_order_dmn,voltron_order(:,3)))
title('Sig. Permuted Coefs FC DMN ordered vs peak drs')
colormap(Color_Salmon_Blue)
hold on
line([47,47],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([117,117],[0,502],'Color','black','LineWidth',1) %somat mot a/b
line([169,169],[0,502],'Color','black','LineWidth',1) %dorsatten
line([220,220],[0,502],'Color','black','LineWidth',1) %sal ventatten
line([244,244],[0,502],'Color','black','LineWidth',1) % limbic
line([305,305],[0,502],'Color','black','LineWidth',1) %all Cont A-C
line([384,384],[0,502],'Color','black','LineWidth',1) %all default nets
line([400,400],[0,502],'Color','black','LineWidth',1) %temp par net
line([414,414],[0,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([430,430],[0,502],'Color','black','LineWidth',1) %thal
line([454,454],[0,502],'Color','black','LineWidth',1) %basal gang.
line([482,482],[0,502],'Color','black','LineWidth',1)
%order DMN into regional grps
%order DMN into regional grps
line([0,502],[9.5,9.5],'Color','black','LineWidth',1) %IPL
line([0,502],[25.5,25.5],'Color','black','LineWidth',1) %PFCd
line([0,502],[37.5,37.5],'Color','black','LineWidth',1) %PFCm
line([0,502],[47.5,47.5],'Color','black','LineWidth',1) %PFCl/v
line([0,502],[59.5,59.5],'Color','black','LineWidth',1) %pCUNPCC
line([0,502],[69.5,69.5],'Color','black','LineWidth',1) %Temp lobe
line([0,502],[74.5,74.5],'Color','black','LineWidth',1) %DefaultC_Rsp
%% Extracted for Chord Plots 
%create variable for chord plot
reorder_fc_dmn = mat_sig_coefs_lm(voltron_order(:,3),new_order_dmn); %ordered schaef x ordered dmn
reorder_fc_dmn = sig_diff_fc_dmn(voltron_order(:,3),new_order_dmn);
writematrix(reorder_fc_dmn,'ordered_fc_dmn_lm.csv')
%grouped into networks for schaef
positive_count = zeros(size(reorder_fc_dmn,2),1);
negative_count = zeros(size(reorder_fc_dmn,2),1);
% Loop through each row to count positive and negative values
for i = 1:size(reorder_fc_dmn, 1)
    positive_count(i) = sum(reorder_fc_dmn(i, :) > 0);
    negative_count(i) = sum(reorder_fc_dmn(i, :) < 0);
end
% take the sum of these within the sub-network of the DMN
remove_zero_fc_dmn = reorder_fc_dmn;
remove_zero_fc_dmn(remove_zero_fc_dmn==0)=NaN;

avg_fc_dmn_nets_order = zeros(35,79);
%avg. across each network for 
for i=1:35
    avg_fc_dmn_nets_order(i,:)=nanmean(remove_zero_fc_dmn((reorder_schaef_nets(:,1)==i),:),1);
end
% avg. across each DMN sub-ROIs
avg_fc_dmn_nets_order2 = zeros(35,8);
for i=1:8
    avg_fc_dmn_nets_order2(:,i)=nanmean(avg_fc_dmn_nets_order(:,(new_order_dmn(:,2)==i)),2);
end

%convert NaN back to 0
avg_fc_dmn_nets_order2(isnan(avg_fc_dmn_nets_order2))=0;
writematrix(avg_fc_dmn_nets_order2,'avg_ordered_fc_dmn_lm.csv')

% create for only 7 networks of brain
load('ordered_schaef_nets_dmn_regions.mat')

remove_zero_fc_dmn_order = mat_sig_coefs_lm(ordered_schaef_dmn(:,2),new_order_dmn(:,2));


remove_zero_fc_dmn(mat_coef==0)=NaN;
for i=1:33
    avg_fc_dmn_nets_order(i,:)=nanmean(remove_zero_fc_dmn((ordered_schaef_dmn(:,2)==i),:),1);
end

% create chord diagram - 
sig_diff_fc_dmn_nonzero = sig_diff_fc_dmn;
sig_diff_fc_dmn_nonzero(sig_diff_fc_dmn==0)=NaN;
avg_across_subnets = zeros(8,502);
for i=1:8
    avg_across_subnets(i,:)= mean(sig_diff_fc_dmn_nonzero(:,order_subnets(:,2)==i),2,'omitmissing');
end
avg_across_subnets(isnan(avg_across_subnets))=0;
avg_across_subnets = avg_across_subnets';
writematrix(avg_across_subnets,'avg_across_dmn_subnets_wholebrain.csv');


%% Create Sub-Network Plots
DMN_subnets = zeros(400,1);
DMN_subnets(:,2)=zeros;

RB_surf_schaef(DMN_subnets(:,2),'')
load('dmn_order_new.mat')
a = DMN_subnets([149:194,358:390],:); %just DMN subnets
[b,order_subnets] = sort(a(:,2)); % core, med temp, dors/med pfc ordering




%% Figures
%load colors
load("colormaps2.mat")
load("schaef_order.mat") %loads in order for voltron into networks
% 1. Permuted significant diff. 
figure
set(gcf,'color','w')
subplot(1,2,1)
%imagesc(sig_fc_del_dmn(voltron_order(:,3),:))
imagesc(sig_fc_del_dmn2(voltron_order(:,3),:))
ylabel('ROIs')
xlabel('DMN Nodes')
title('Sig. Avg. DMN FC Delirium')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum
subplot(1,2,2)
%imagesc(sig_fc_health_dmn(voltron_order(:,3),:))
imagesc(sig_fc_health_dmn2(voltron_order(:,3),:))
ylabel('ROIs')
xlabel('DMN nodes')
title('Sig. Avg. DMN FC Health')
colormap(Color_Salmon_Blue)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum

figure
set(gcf,'color','w')
imagesc(sig_diff_fc_dmn(voltron_order(:,3),:))
ylabel('ROIs')
xlabel('DMN nodes')
title('Sig. Avg. Diff DMN FC Del - Health')
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %


figure
set(gcf,'color','w')
%imagesc(reshape_sig_rho_FC_DMN(voltron_order(:,3),:))
imagesc(reshape_sig_permute_corr(voltron_order(:,3),:))
ylabel('ROIs')
title('Sig. Permuted Rho FC DMN vs peak drs')
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)


figure
set(gcf,'color','w')
subplot(1,2,1)
imagesc(sig_fc_del_dmn([149:194,358:390],:))
ylabel('ROIs')
xlabel('DMN Nodes')
title('Sig. DMN FC Delirium')
hold on
line([0,79],[45,45],'Color','black','LineWidth',1) %hemisphere line
line([45,45],[79,0],'Color','black','LineWidth',1) %hemisphere line
subplot(1,2,2)
imagesc(sig_fc_health_dmn([149:194,358:390],:))
ylabel('ROIs')
xlabel('DMN nodes')
title('Sig. DMN FC Healthy')
colormap(Color_Salmon_Blue)
hold on
line([0,79],[45,45],'Color','black','LineWidth',1) %hemisphere line
line([45,45],[79,0],'Color','black','LineWidth',1) %hemisphere line




%% Secondary Analysis - DMN avg time-series, then correlated
cd 'C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\timeseries\schaef_400'
filename=dir('*_timeseries.mat');
fc_all_avg_dmn = zeros(127,425,425);
%DMN ids
[149:194,358:390]
for ii=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(ii).name; %filename of time-series data
    split = strsplit(subject_file,'_');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum)
    ses = split(2); %session ses-..
    ses = cell2mat(ses);
    load([subject_file]); %load in time-series; ts variable
    ts = ts(6:end -5,:); %remove first/last 5 time-points for noise
    %functional connectivity just to avg. time-series of L & R Hemispheres
    avg_dmn_L = mean(ts(:,149:194),2);
    avg_dmn_R = mean(ts(:,358:390),2);
    ts_remove_dmn = ts;
    ts_remove_dmn(:,[149:194,358:390]) =[];
    ts2 = [ts_remove_dmn,avg_dmn_L,avg_dmn_R];
    ts_corr = corr(ts2); %ROI x ROI matrix
    fc_all_avg_dmn(ii,:,:)=ts_corr; %all subjects FC into one
end




load('C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Demographic_Data\V1_MRI_all_demos.mat');
sub_remove = find(subjectV1MRIdata.overall_delirious_ever=='NA');
subjectV1MRIdata(sub_remove,:)=[];
bin_delirium_all=table2array(subjectV1MRIdata(1:size(subjectV1MRIdata,1),"bin_delirium"));
delirium_sub = find(bin_delirium_all==1);
health_sub =bin_delirium_all~=1; 
health_sub = find(health_sub==1);
fc_all_avg_dmn(sub_remove,:,:)=[];

%remove subjects without peak del. data
fc_del_dmn = fc_all_avg_dmn(delirium_sub,:,:);
fc_health_dmn = fc_all_avg_dmn(health_sub,:,:);

% group-wise permutation significant of PC edges
nROIs = size(fc_del_dmn_flat,2);
fc_del_dmn_flat = reshape(fc_del_dmn,31,425*425);
fc_health_dmn_flat = reshape(fc_health_dmn,89,425*425);

for i=1:nROIs
   [sig_fc_avg_dmn(i,:),pval_fc_avg_dmn(i,:)]=perm_code(fc_del_dmn_flat(:,i),fc_health_dmn_flat(:,i),1000);
end


sig_fc_avg_dmn_del = reshape(sig_fc_avg_dmn,425,425).*squeeze(mean(fc_del_dmn));
sig_fc_avg_dmn_health = reshape(sig_fc_avg_dmn,425,425).*squeeze(mean(fc_health_dmn));

sig_diff_fc_avg_dmn = sig_fc_avg_dmn_del - sig_fc_avg_dmn_health;

%make new voltron order
order_no_dmn = voltron_order(:,2);
order_no_dmn([149:194,358:390])=[];
order_no_dmn([424,425],1)=36;
[sort,ordered_dmn] = sort(order_no_dmn); %reorders without dmn network included



figure
set(gcf,'color','w')
subplot(1,2,1)
%imagesc(sig_fc_del_dmn(voltron_order(:,3),:))
imagesc(sig_fc_avg_dmn_del(ordered_dmn,:))
ylabel('ROIs')
xlabel('DMN Nodes')
title('Sig. Avg. DMN FC Delirium')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum
subplot(1,2,2)
%imagesc(sig_fc_health_dmn(voltron_order(:,3),:))
imagesc(sig_fc_avg_dmn_health(ordered_dmn,:))
ylabel('ROIs')
xlabel('DMN nodes')
title('Sig. Avg. DMN FC Health')
colormap(Color_Salmon_Blue)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum

figure
set(gcf,'color','w')
imagesc(sig_fc_avg_dmn_health(voltron_order(:,3),:))
ylabel('ROIs')
xlabel('DMN nodes')
title('Sig. Avg. Diff DMN FC Del - Health')
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %
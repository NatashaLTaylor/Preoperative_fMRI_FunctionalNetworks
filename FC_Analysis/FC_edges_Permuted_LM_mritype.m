%% Update FC relationship with delirium severity (peak overall DRS) control for MRI type

%load in the data
cd('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/FC/lmmodel/mri_type_covar/')
%load data local
cd('/Users/ntaylor/Desktop/mri_type_covar/')

data = zeros(125751,2);
%loading in data etc..
for dd=1:126
    filename = sprintf("permute_lm_coef_pval_interact_arrayjob%d.csv",dd);
    % load in the file
    permute_coefs = readtable(filename);
    %change format of coef into values only
    column_data1 = permute_coefs{:,1}; %coefs
    column_data2 = permute_coefs.p_value; %pvalues
    %now need to do this for both columns
    new_coef = zeros(size(permute_coefs,1),2);
    new_coef(:,1)= column_data2;
    new_coef(:,2)= column_data1; %converted coef fc only
%load into one large matrix
    if dd ~=126
        strIdx = (dd-1)*1000+1;
        data(strIdx:strIdx+999,:)=new_coef;
    else %126 less than 1000 places
        strIdx = (dd-1)*1000+1; %this is for 126 
        data(strIdx:strIdx+750,:)=new_coef;
    end
    sprintf("completed arrayjob no. %d",dd)
end

data = zeros(125751,2);% Preallocate data matrix

%For Files that are missing
for dd = 1:126
    filename = sprintf("permute_lm_coef_pval_interact_arrayjob%d.csv", dd);
    
    if exist(filename, 'file') % Check if the file exists
        % Load in the file
        permute_coefs = readtable(filename);
        
        % Change format of coef into values only
        column_data1 = permute_coefs{:,1}; % Coefs
        column_data2 = permute_coefs.p_value; % P-values
        
        % Initialize new_coef with the correct size
        new_coef = zeros(size(permute_coefs, 1), 2);
        new_coef(:, 1) = column_data2;
        new_coef(:, 2) = column_data1; % Converted coef fc only
    else
        % File does not exist, set new_coef to zeros
        if dd ~= 126
            new_coef = zeros(1000, 2); % 1000 rows for most cases
        else
            new_coef = zeros(751, 2); % 751 rows for the last case
        end
    end
    
    % Load into one large matrix
    if dd ~= 126
        strIdx = (dd-1)*1000 + 1;
        data(strIdx:strIdx+999, :) = new_coef;
    else
        strIdx = (dd-1)*1000 + 1; % This is for 126
        data(strIdx:strIdx+750, :) = new_coef;
    end
    
    fprintf("Completed arrayjob no. %d\n", dd);
end


% column 1 p-value fc only, coef fc only, column 3 p=value interaction,
% column 4 coef interaction

% look at only the significant values

% calculate for control for covariate
sig_permute_fc = double(data(:,1)<0.05);
sig_coef_fc = sig_permute_fc.*data(:,2); %multiply the beta coefs with the p-value
%sum of all significant edges
sum(sig_permute_fc)
mat_sig_coef_fc = matify(sig_coef_fc,502);

mat_coef_fc = matify(data(:,2),502);

%% Figures
%plot these

load('schaef_order.mat')
load('colormaps2.mat')

figure
set(gcf,'color','w')
%subplot(1,2,1)
%imagesc(mat_coef_fc(voltron_order(1:400,3),voltron_order(1:400,3)))
imagesc(mat_sig_coef_fc(voltron_order(1:400,3),voltron_order(1:400,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Peak DRS ~ FC edges + MRI Scanner')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,400],'Color','black','LineWidth',1)
line([1,400],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([1,400],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([1,400],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([1,400],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([1,400],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([1,400],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue2)
%remove ticks
xticks([])
yticks([])
xticklabels([])
yticklabels([])


%plot just fc only relationship
%figure
%set(gcf,'color','w')
subplot(1,2,2)
imagesc(mat_sig_coef_fconly(voltron_order(1:400,3),voltron_order(1:400,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Peak DRS ~ FC edges only')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,400],'Color','black','LineWidth',1)
line([1,400],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([1,400],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([1,400],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([1,400],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([1,400],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([1,400],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue)
%remove ticks
xticks([])
yticks([])
xticklabels([])
yticklabels([])

%% Plot half of the edges + the average relationship
%just plot the lower triangle 
order_mat_sig_coef_fc = mat_sig_coef_fc(voltron_order(:,3),voltron_order(:,3));

tril = tril(order_mat_sig_coef_fc);

figure
set(gcf,'Color','w')
imagesc(tril(1:400,1:400))
hold on
line([1,400],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,400],'Color','black','LineWidth',1)
line([1,400],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([1,400],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([1,400],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([1,400],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([1,400],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([1,400],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue)
xticklabels([])
yticklabels([])
% calculate avg. across








%% Plot of alll 
figure
set(gcf,'color','w')
%subplot(1,2,1)
imagesc(mat_sig_coef_interact(voltron_order(:,3),voltron_order(:,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Cortical ROIs')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
%subcortical
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
%rest to 502 are asending arousal
colormap(Color_Salmon_Blue)










    %calculated for just relevant brain relationship
    sig_permute = double(permute_coefs(:,1)<0.05);
    sig_coef = sig_permute.*permute_coefs(:,2); %multiply the beta coefs with the p-value
    length_coef = length(sig_coef);
    sig_coef(length_coef,:)=[]; %last variable is the covariate coef

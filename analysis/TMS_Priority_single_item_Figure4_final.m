% TMS_Priority_single_item_wmChooseTopoPriority.m
% Compare two-item high priority data from this study to single-item
% data from indepenent data sets (partial subject overlap)
% Produces Figure 4.
%
% Hallenbeck, Tardiff, et al., JNeurosci, 2025

clear all;close all

%% setup
%matlab defaults for appearance
set(groot,'DefaultFigureColor',[1 1 1])

set(groot,'DefaultAxesFontName','Arial')
set(groot,'DefaultLegendFontName','Arial');
set(groot,'DefaultTextFontName','Arial');
set(groot,'DefaultAxesFontSize',7); %15)

set(groot,'DefaultLineLineWidth',1.75)
set(groot,'DefaultErrorBarLineWidth',2.5)
set(groot,'DefaultErrorBarMarkerSize',12)
set(groot,'DefaultErrorBarMarkerFaceColor','auto')
set(groot,'DefaultErrorBarCapSize',0)
set(groot,'DefaultScatterLineWidth',1.75)
set(groot,'DefaultScatterMarkerFaceColor',[1 1 1]);
set(groot,'DefaultLegendBox','off')


root = '../';
%addpath(genpath(fullfile(root,'util')));
data_dir = fullfile(root,'data');
fig_dir = fullfile(root,'figs_raw');

subj_valid_priority = [1 2 3 4 5 6 7 11 12 16 17 21 22 24];
hemi = {'LVF','RVF'};


%data files
two_data_file = fullfile(data_dir,'TMS_priority_MGS_data_11-Oct-2024_allconds.mat');
single_data_file = fullfile(data_dir,'TMS_singleitem_MGS_data_07-Mar-2024.mat');
mrugank_data_file = fullfile(data_dir,'mrugank_tms_summary_byhemi.csv');

%read in data and harmonize variable names and select only noTMS
%conditions,only high priority for two-item

%two-item
F = load(two_data_file);
data_2item = F.mgs_data;

%make sure subjects excluded
data_2item = data_2item(ismember(data_2item.subject,subj_valid_priority),:);
assert(all(ismember(data_2item.subject,subj_valid_priority)),'Incorrect subjects detected!');

data_2item = data_2item(strcmp(data_2item.condition,'noTMS') & ...
    strcmp(data_2item.priority,'high'),:); 
data_2item.study = cellstr(repmat('two_item',height(data_2item),1));


%single-item Topo/MGSMap
F = load(single_data_file);
data_1item = F.mgs_data;
clear F
data_1item.study = cellstr(repmat('one_item',height(data_1item),1));

%single-item Mru
data_1itemM = readtable(mrugank_data_file);
data_1itemM = data_1itemM(strcmp(data_1itemM.TMS_time,'notms'),:);
data_1itemM.study = cellstr(repmat('one_itemM',height(data_1itemM),1));
data_1itemM.hemi_raw = double(...
    (strcmp(data_1itemM.hemistimulated_first,'Left') & data_1itemM.instimVF) | ...
    (strcmp(data_1itemM.hemistimulated_first,'Right') & ~data_1itemM.instimVF)) + 1;
data_1itemM.hemi = hemi(data_1itemM.hemi_raw)';
data_1itemM = renamevars(data_1itemM,{'ierr_mean'},{'mean_i_sacc_err'});

%% combine data, create means
data_all = [data_2item(:,{'subject','study','hemi','i_sacc_err'});...
    data_1item(:,{'subject','study','hemi','i_sacc_err'})];
data_all_ave = varfun(@mean,data_all,'InputVariables',{'i_sacc_err'},...
    'GroupingVariables',{'subject','study','hemi'});
data_all_ave.GroupCount = [];
data_all_ave = [data_all_ave; data_1itemM(:,{'subject','study','hemi','mean_i_sacc_err'})];
data_all_ave.load = strcmp(data_all_ave.study,'two_item') + 1;

%% spot-checks
%find common subjects between two- and one-item data
subj_21 = intersect(data_2item.subject,data_1item.subject); %in two and one
subj_21M = intersect(data_2item.subject,data_1itemM.subject); %in two and oneM
subj_211M = intersect(subj_21,subj_21M); %in two, one, and oneM
subj_11M = intersect(data_1item.subject,data_1itemM.subject); %in one and oneM


%% plot averages


%create average for between-subjects version
%average over subjects w/ multiple data points per load
data_all_ave2 = varfun(@mean,data_all_ave,"InputVariables",'mean_i_sacc_err',...
    'GroupingVariables',{'subject','load','hemi'});
data_all_ave2 = renamevars(data_all_ave2,'mean_mean_i_sacc_err','mean_i_sacc_err');
%group averages
data_all_ave2_ave = summary_stats(data_all_ave2,'mean_i_sacc_err',{'load','hemi'});

%let's refind the subjects w/ both one and two item data here...was janky
%the way I was doing it before
subj_onetwo = intersect(data_all_ave2.subject(data_all_ave2.load==1), ...
    data_all_ave2.subject(data_all_ave2.load==2));

data_all_12 = data_all_ave2(ismember(data_all_ave2.subject,subj_onetwo),:);
data_all_12_ave = summary_stats(data_all_12,'mean_i_sacc_err',{'hemi','load'});

%we're just going to plot rh
plot_hemi = {'RVF'}; %,'LVF'};

subplot = @(m,n,p) subtightplot(m,n,p,.12,[0.11,0.05],[.09,0.015]);
for h=1:length(plot_hemi)
    f = figure('units','centimeters','Position',[20 20 16 12]);

    %first plot all subjects 
    this_data_ave = data_all_ave2(strcmp(data_all_ave2.hemi,plot_hemi{h}),:);
    this_data_ave_all = data_all_ave2_ave(strcmp(data_all_ave2_ave.hemi,plot_hemi{h}),:);

    y_max = max(this_data_ave.mean_i_sacc_err)*1.05;
        
    subplot(1,2,1);hold on
    %rng(h+1)
    bar(this_data_ave_all.load,this_data_ave_all.mean_mean_i_sacc_err,'FaceColor',[0.5 0.5 0.5],...
        'EdgeAlpha',1,'BarWidth',0.7);
    %swarmchart(this_data_ave.load,this_data_ave.mean_i_sacc_err,...
    scatter(this_data_ave.load,this_data_ave.mean_i_sacc_err,...
        'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.5,...
        'MarkerFaceAlpha',0.5)
    errorbar(this_data_ave_all.load,this_data_ave_all.mean_mean_i_sacc_err,...
        this_data_ave_all.sem_mean_i_sacc_err,'.k');%'ok')
    xticks([1 2]);
    xlabel('load')
    %xlim([0.25 2.75])
    xlim([0.5 2.5])
    ylim([0 y_max])
    ylabel('Error (\circ)')
    

    %now plot within-subject comparison (actually plotting this first but
    %didn't rearrange)
    this_data_12 = data_all_12(strcmp(data_all_12.hemi,plot_hemi{h}),:);
    this_data_12_ave = data_all_12_ave(strcmp(data_all_12_ave.hemi,plot_hemi{h}),:);

    this_data_12W = unstack(this_data_12,'mean_i_sacc_err','load',...
        'GroupingVariables','subject','NewDataVariableNames',{'load1','load2'});
    load_mat = repmat([1 2],height(this_data_12W),1);
    err_mat = this_data_12W{:,{'load1','load2'}};

    subplot(1,2,2);hold on;
    
    bar(this_data_12_ave.load,this_data_12_ave.mean_mean_i_sacc_err,'FaceColor','r',...
        'EdgeAlpha',1,'BarWidth',0.7);
    plot(load_mat',err_mat','-','Color',[0 0 0 0.5],'LineWidth',1.5)
    errorbar(this_data_12_ave.load,this_data_12_ave.mean_mean_i_sacc_err,...
        this_data_12_ave.sem_mean_i_sacc_err,'.k');%'ok')
    xlim([0.5 2.5])
    xticks([1 2]);
    xlabel('load')
    ylim([0 y_max])
    %ylabel('Error (\circ)')

    %save figure
    exportgraphics(f,fullfile(fig_dir,sprintf('singleitem_comp_%s_%s.pdf', ...
        plot_hemi{h},date())), ...
        'ContentType','vector')

end


%% stats

%overall difference between one and two item?
data_all_ave2_RVF = data_all_ave2(strcmp(data_all_ave2.hemi,'RVF'),:); 
[t,stats,orig,param] = ...
    permtest_t2(data_all_ave2_RVF.mean_i_sacc_err(data_all_ave2_RVF.load==1,:),...
    data_all_ave2_RVF.mean_i_sacc_err(data_all_ave2_RVF.load==2,:))

Mdiff = data_all_ave2_ave.mean_mean_i_sacc_err(strcmp(data_all_ave2_ave.hemi,'RVF') & ...
    data_all_ave2_ave.load==1,:) - ...
    data_all_ave2_ave.mean_mean_i_sacc_err(strcmp(data_all_ave2_ave.hemi,'RVF') & ...
    data_all_ave2_ave.load==2,:)


%within subject comparison for those in both datasets
data_all_12W_RVF =  unstack(data_all_12(strcmp(data_all_12.hemi,'RVF'),:),...
    'mean_i_sacc_err','load',...
    'GroupingVariables','subject','NewDataVariableNames',{'load1','load2'});

[t,stats,orig,param] = ...
    permtest_t(data_all_12W_RVF.load1,data_all_12W_RVF.load2)

Mdiff_within = mean(data_all_12W_RVF.load1) - ...
   mean(data_all_12W_RVF.load2)





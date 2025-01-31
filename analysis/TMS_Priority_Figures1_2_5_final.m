% Figure 1 example scatter and noTMS effects over all data plots/analysis
% Fig 2 behavioral effect of TMS analysis PCS vs. noTMS
% Fig 5 behavioral effect of TMS analysis IPS2 vs. noTMS
% Fig 2 PCS efield strength x hemi
% Contains graphing code and post-hoc analyses
% For ANOVA see mgs_analysis.R
%
% Hallenbeck, Tardiff, et al., JNeurosci, 2025


close all;clear all

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


addpath(genpath('../util'));
data_dir = '../data'; 
fig_dir = '../figs_raw';
data_file = 'TMS_priority_MGS_data_11-Oct-2024_allconds.mat';

% FOR INITIAL SACCADE ANALYSIS % mae %%%%%%%% EXCLU 8,9,20
subj_valid = [1 2 3 4 5 6 7 11 12 16 17 21 22 24];
subj_excl = [8 9 20];
subj = sort([subj_valid subj_excl]);

TMS_cond = {'noTMS','l_spcs','l_ips2'}; 
TMS_cond_names = {'none','lSPCS','lIPS2'}; 
pcond = {'high','low'};
cond_lines = {'-','--',':'};

%load already preprocessed data
F = load(fullfile(data_dir,data_file));
all_data = F.mgs_data;
clear F

%we only want noTMS and l_spcs in these analyses:
all_data = all_data(ismember(all_data.condition,TMS_cond),:); %we only want noTMS and IPS here

%check that data matches what we expect
assert(isequal(subj,unique(all_data.subject)'),'Subject list does not match data!')
assert(isequal(sort(TMS_cond),unique(all_data.condition)'),'Condition list does not match data!')
assert(isequal(sort(pcond),unique(all_data.priority)'),'Condition list does not match data!')
assert(~any(isnan(all_data.i_sacc_err)),'Missing data detected!')

%for figures
mgs_meas = {'mean_i_sacc_err','std_i_sacc_err','mean_i_sacc_rt'};
mgs_ylabel = {'Error (\circ)','Errror SD (\circ)','RT (s)'};
plot_ylabel = struct();
for m=1:length(mgs_meas)
    plot_ylabel.(mgs_meas{m}) = mgs_ylabel{m};
end
cond_color = lines(2);

%% compute averages for plotting/analysis

%compute averages
mgs_vars = {'i_sacc_err'}; 
mgs_rtvars = {'i_sacc_rt'}; 
mgs_meanvars = [strcat('mean_',[mgs_vars mgs_rtvars])]; 

%we compute mean error and variance of error as complementary
data_ave_full = summary_stats(all_data,[mgs_vars mgs_rtvars], ...
    {'subject','condition','priority','hemi'});
data_ave_full_all = summary_stats(data_ave_full,...
    mgs_meanvars,{'condition','priority','hemi'});

%also compute the prioritization effect (low - high) for visualization
data_aveW_full = unstack(data_ave_full(:,[{'subject','condition','priority','hemi'},mgs_meanvars]),...
    mgs_meanvars,'priority', ...
    'GroupingVariables',{'subject','condition','hemi'});
peffect_inputs = {'mean_i_sacc_err','mean_i_sacc_rt'};
peffect_vars = strcat(peffect_inputs,'_diff');

for i=1:length(peffect_vars)
    data_aveW_full.(peffect_vars{i}) = ...
        data_aveW_full.([peffect_inputs{i} '_low']) - data_aveW_full.([peffect_inputs{i} '_high']);
end

%% noTMS analysis over all subjects

data_ave_full_noTMS = data_ave_full(strcmp(data_ave_full.condition,'noTMS'),:);
data_ave_full_noTMS_all = data_ave_full_all(strcmp(data_ave_full_all.condition,'noTMS'),:);

plot_meas = {'mean_i_sacc_err','mean_i_sacc_rt'};
plot_hemi = {'RVF','LVF'}; %{'RVF','LVF'};
%num_ticks = 4;

subplot = @(m,n,p) subtightplot(m,n,p,.08,[0.11,0.05],[.08,0.025]);
for h=1:length(plot_hemi)
    f=figure('units','centimeters','Position',[20 20 16 12]);
    for m=1:length(plot_meas)
        subplot(1,length(plot_meas),m);hold on;
        
        %plot high vs. low priority
        this_pdata = data_ave_full_noTMS(strcmp(data_ave_full_noTMS.hemi,plot_hemi{h}),:);
        this_pdata.priority = categorical(this_pdata.priority,...
            pcond);
        this_subj = unique(this_pdata.subject);
        this_pcond = unique(this_pdata.priority);
        this_ymax = round(max(this_pdata.(plot_meas{m})) * 1.05,1);
        this_ylim = [0 this_ymax];

        this_pdata_ave = data_ave_full_noTMS_all(strcmp(data_ave_full_noTMS_all.hemi,...
            plot_hemi{h}),:);
        this_pdata_ave.priority = categorical(this_pdata_ave.priority,...
            pcond);

        %yline(0,'--k')
        %scatter(double(this_pdata.priority),this_pdata.(plot_meas{m}),65, ...
        %    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
        for ss=1:length(this_subj)
            this_pdata_subj = this_pdata(this_pdata.subject==this_subj(ss),:);
            plot(double(this_pdata_subj.priority),this_pdata_subj.(plot_meas{m}), ...
                '-','Color',[0 0 0 0.5],'LineWidth',1);
        end
        errorbar(double(this_pdata_ave.priority),...
            this_pdata_ave.(['mean_' plot_meas{m}]),...`
            this_pdata_ave.(['sem_' plot_meas{m}]),'-ok');
        
        xticks([1 2])
        xticklabels(this_pcond)
        xlim([0.5 2.5])
        ylim(this_ylim)
        %yticks(this_ylim(1):(round(this_ylim(2) - this_ylim(1),1)/num_ticks):this_ylim(2))
        ylabel(plot_ylabel.(plot_meas{m}))
    end

    %save figure
    exportgraphics(f,fullfile(fig_dir,sprintf('noTMS_behave_%s_%s.pdf', ...
        plot_hemi{h},date())), ...
        'ContentType','vector')
    
end


%% noTMS analysis
data_ave_full_noTMSW_RVF = ...
    data_aveW_full(strcmp(data_aveW_full.hemi,'RVF') &...
    strcmp(data_aveW_full.condition,'noTMS'),:);

data_ave_full_noTMSW_LVF = ...
    data_aveW_full(strcmp(data_aveW_full.hemi,'LVF') &...
    strcmp(data_aveW_full.condition,'noTMS'),:);


%RVF acc and RT
[t,stats,orig,param] = permtest_t(data_ave_full_noTMSW_RVF.mean_i_sacc_err_high,...
    data_ave_full_noTMSW_RVF.mean_i_sacc_err_low)


[t,stats,orig,param] = permtest_t(data_ave_full_noTMSW_RVF.mean_i_sacc_rt_high,...
    data_ave_full_noTMSW_RVF.mean_i_sacc_rt_low)


%LVF acc and RT (supp)
[t,stats,orig,param] = permtest_t(data_ave_full_noTMSW_LVF.mean_i_sacc_err_high,...
    data_ave_full_noTMSW_LVF.mean_i_sacc_err_low)


[t,stats,orig,param] = permtest_t(data_ave_full_noTMSW_LVF.mean_i_sacc_rt_high,...
    data_ave_full_noTMSW_LVF.mean_i_sacc_rt_low)


%determine who doesn't show priority effect (excluding based on RVF only!)
nop_subj = data_ave_full_noTMSW_RVF(data_ave_full_noTMSW_RVF.mean_i_sacc_err_diff <= 0,:)
nop_subj_LVF = data_ave_full_noTMSW_RVF(data_ave_full_noTMSW_LVF.mean_i_sacc_err_diff <= 0,:)

assert(isequal(nop_subj.subject,subj_excl'))

%% exclude subjects who don't show priority effect for TMS analysis

%remove subjects who don't show the effect from data structures 
data_ave = data_ave_full(~ismember(data_ave_full.subject,nop_subj.subject),:);
data_aveW = data_aveW_full(~ismember(data_aveW_full.subject,nop_subj.subject),:);

%compute a new data ave and data aveW and all data (or just subset all data as needed)
data_ave_all = summary_stats(data_ave,...
    mgs_meanvars,{'condition','priority','hemi'});
data_aveW_all = summary_stats(data_aveW,peffect_vars,{'condition','hemi'});


%% plot saccade endpoints for example subject
data_aveW_noTMS_RVF = data_aveW(strcmp(data_aveW.hemi,'RVF') & ...
    strcmp(data_aveW.condition,'noTMS'),:);

%find subject with largest noTMS priority effect
[~,max_pr_I] = max(data_aveW_noTMS_RVF.mean_i_sacc_err_diff);
ex_subj = data_aveW_noTMS_RVF.subject(max_pr_I);

%figure out limits
targ_ecc = 9;
ex_subj_data = all_data(strcmp(all_data.hemi,'RVF') & ...
        strcmp(all_data.condition,'noTMS') & all_data.subject==ex_subj,:);
ex_min_raw = min(ex_subj_data.i_sacc) - [targ_ecc 0];
ex_max_raw = max(ex_subj_data.i_sacc) - [targ_ecc 0];
ex_min = floor(ex_min_raw);
ex_max = ceil(ex_max_raw);

x_range = [-7 7]; %[ex_min(1) ex_max(1)] + [0 1]; %[-1 17.5] - targ_ecc;
y_range = [-6 6]; %[ex_min(2) ex_max(2)] + [0 1]; %[-8,8];

subplot = @(m,n,p) subtightplot(m,n,p,[0.01,.02],[.15 .1],.1);
f = figure('Units','centimeters','Position',[0 0 30 15]);
for i=1:length(pcond)
    
    this_subj = all_data(strcmp(all_data.hemi,'RVF') & strcmp(all_data.priority,pcond(i)) &...
        strcmp(all_data.condition,'noTMS') & all_data.subject==ex_subj,:);
    
    this_subj.i_sacc(:,1) = this_subj.i_sacc(:,1) - targ_ecc; %target centered at 9 ecc
    
    subplot(1,2,i);hold on;

    xline(0,'--','Alpha',0.6,'LineWidth',1)
    yline(0,'--','Alpha',0.6,'LineWidth',1)
    
    
    xgrid=linspace(x_range(1),x_range(2),100);
    ygrid=linspace(y_range(1),y_range(2),100);
    [x1,y1] = meshgrid(xgrid, ygrid);
    % Perform kernel density estimate
    % [x y] is actual data, xi is the desired grid points to evaluate
    % fd is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
    xi = [x1(:) y1(:)];
    [fd,ep]=ksdensity(this_subj.i_sacc,xi,"Bandwidth",1.5); %,1.5); % remove the outputs to see a 3D plot of the distribution
    % format data in matrix for contourf and plot
    X = reshape(ep(:,1),length(xgrid),length(ygrid));
    Y = reshape(ep(:,2),length(xgrid),length(ygrid));
    Z = reshape(fd,length(xgrid),length(ygrid));
    
    %contour(X,Y,Z,6,'LineWidth',1)
    %[M,c] = contourf(X,Y,Z,6,'LineWidth',1);
    [M,h]=contourf(X,Y,Z,6,'LineStyle','none');
    %h.FaceAlpha = 0.5;  
    colormap(brewermap([],"YlOrRd"))
    
    
    %c.ZData(c.ZData < c.LevelList(1)) = NaN;
    %ksdensity(this_subj.i_sacc,'PlotFcn','contour','Bandwidth',1.5);

    %plot(this_subj.i_sacc(:,1),this_subj.i_sacc(:,2),'k.');
    scatter(this_subj.i_sacc(:,1),this_subj.i_sacc(:,2),30,...
        'MarkerFaceColor','white','MarkerFaceAlpha',1.0,...
        'MarkerEdgeColor','black','MarkerEdgeAlpha',1.0)

    axis equal
    box off

    yticks([-4 0 4])
    xticks([-8 -4 0 4])

    if i > 1
        yticklabels([]);
    end

end

exportgraphics(f,fullfile(fig_dir,sprintf('sacc_scatter_%s.pdf',date())), ...
    'ContentType','vector')


%% Plot MGS error x condition (notms,spcs)
plot_meas = {'mean_i_sacc_err','mean_i_sacc_rt'}; 
plot_hemi = {'RVF','LVF'}; 
TMS_cond_plot2 = 1:2;

subplot = @(m,n,p) subtightplot(m,n,p,.09,[0.11,0.05],[.09,0.015]);
%clear subplot
for h=1:length(plot_hemi)
    
    %pull out hemi data
    this_data_hemi = data_ave_all(strcmp(data_ave_all.hemi,plot_hemi{h}),:);
    this_data_subj_hemi = data_ave(strcmp(data_ave.hemi,plot_hemi{h}),:);
    
    this_subj = unique(this_data_subj_hemi.subject);
    
    %get per hemi limits
    meas_maxes = struct();
    for m=1:length(plot_meas)
        this_meas = plot_meas{m};
        this_fudge = (max(this_data_subj_hemi.(this_meas)) - min(this_data_subj_hemi.(this_meas)))*0.05;
        meas_maxes.(plot_meas{m}) = ...
            max(this_data_subj_hemi.(this_meas))+this_fudge;
    end

    f=figure('units','centimeters','Position',[20 20 16 12]);
    for m=1:length(plot_meas)

        subplot(1,length(plot_meas),m);hold on;
        for c=TMS_cond_plot2 %length(TMS_cond)
            %first plot effect in each condition
            this_data = this_data_hemi(strcmp(this_data_hemi.condition,TMS_cond{c}),:);
            this_data.priority = categorical(this_data.priority);
            errorbar(this_data.priority,this_data.(['mean_' plot_meas{m}]),...
                this_data.(['sem_' plot_meas{m}]),[cond_lines{c} 'ko']); %'-o');

            %now plot individual subject lines
            this_data_subjc = this_data_subj_hemi(strcmp(this_data_subj_hemi.condition,TMS_cond{c}),:);
            this_data_subjc.priority = categorical(this_data_subjc.priority);
            for ss=1:length(this_subj)
                this_pdata_subj = this_data_subjc(this_data_subjc.subject==this_subj(ss),:);
                plot(this_pdata_subj.priority,this_pdata_subj.(plot_meas{m}), ...
                    [cond_lines{c} 'k'],'Color',[0 0 0 0.5],'LineWidth',1); 
                    %'-','Color',[cond_color(c,:) 0.5],'LineWidth',1);
            end
        end

        ylabel(plot_ylabel.(plot_meas{m}));
        ylim([0 meas_maxes.(plot_meas{m})]);



    end
    %save figure
    exportgraphics(f,fullfile(fig_dir,sprintf('%s_%s_%s.pdf', ...
        strjoin(plot_meas,'-'),plot_hemi{h},date())), ...
        'ContentType','vector')
    
end


%% Plot MGS error x condition (notms,ips2)
plot_meas = {'mean_i_sacc_err','mean_i_sacc_rt'}; 
plot_hemi = {'RVF','LVF'}; 
TMS_cond_plot3 = [1,3];


subplot = @(m,n,p) subtightplot(m,n,p,.09,[0.11,0.05],[.09,0.015]);
%clear subplot
for h=1:length(plot_hemi)
    
    %pull out hemi data
    this_data_hemi = data_ave_all(strcmp(data_ave_all.hemi,plot_hemi{h}),:);
    this_data_subj_hemi = data_ave(strcmp(data_ave.hemi,plot_hemi{h}),:);
    
    this_subj = unique(this_data_subj_hemi.subject);
    
    %get per hemi limits
    meas_maxes = struct();
    for m=1:length(plot_meas)
        this_meas = plot_meas{m};
        this_fudge = (max(this_data_subj_hemi.(this_meas)) - min(this_data_subj_hemi.(this_meas)))*0.05;
        meas_maxes.(plot_meas{m}) = ...
            max(this_data_subj_hemi.(this_meas))+this_fudge;
    end

    f=figure('units','centimeters','Position',[20 20 16 12]);
    for m=1:length(plot_meas)

        subplot(1,length(plot_meas),m);hold on;
        for c=TMS_cond_plot3 %length(TMS_cond)
            %first plot effect in each condition
            this_data = this_data_hemi(strcmp(this_data_hemi.condition,TMS_cond{c}),:);
            this_data.priority = categorical(this_data.priority);
            errorbar(this_data.priority,this_data.(['mean_' plot_meas{m}]),...
                this_data.(['sem_' plot_meas{m}]),[cond_lines{c} 'ko']); %'-o');

            %now plot individual subject lines
            this_data_subjc = this_data_subj_hemi(strcmp(this_data_subj_hemi.condition,TMS_cond{c}),:);
            this_data_subjc.priority = categorical(this_data_subjc.priority);
            for ss=1:length(this_subj)
                this_pdata_subj = this_data_subjc(this_data_subjc.subject==this_subj(ss),:);
                plot(this_pdata_subj.priority,this_pdata_subj.(plot_meas{m}), ...
                    [cond_lines{c} 'k'],'Color',[0 0 0 0.5],'LineWidth',1); 
                    %'-','Color',[cond_color(c,:) 0.5],'LineWidth',1);
            end
        end

        ylabel(plot_ylabel.(plot_meas{m}));
        ylim([0 meas_maxes.(plot_meas{m})]);



    end
    %save figure
    exportgraphics(f,fullfile(fig_dir,sprintf('%s_%s_%s_ips.pdf', ...
        strjoin(plot_meas,'-'),plot_hemi{h},date())), ...
        'ContentType','vector')
    
end


%% statistical testing of TMS effects
% NOTE: see: MGS_analysis.R for the ANOVA.

%post-hoc paired-sample permutation t-tests 

%easiest to have all conditions as columns for this...
data_noTMS_RVF = data_aveW(strcmp(data_aveW.condition,'noTMS') & ...
    strcmp(data_aveW.hemi,'RVF'),:);
data_lspcs_RVF = data_aveW(strcmp(data_aveW.condition,'l_spcs') & ...
    strcmp(data_aveW.hemi,'RVF'),:);

%data for across-condition comparisons within priority
data_noTMS_RVF_mat = data_noTMS_RVF{:,{'mean_i_sacc_err_high','mean_i_sacc_err_low'}};
data_lspcs_RVF_mat = data_lspcs_RVF{:,{'mean_i_sacc_err_high','mean_i_sacc_err_low'}};

%data for across-priority comparisons within conditon
%notms vs lspcs
data_conds_noTMS_lspcs_high = [data_noTMS_RVF.mean_i_sacc_err_high data_lspcs_RVF.mean_i_sacc_err_high];
data_conds_noTMS_lspcs_low = [data_noTMS_RVF.mean_i_sacc_err_low data_lspcs_RVF.mean_i_sacc_err_low];


% how many subjects show greater error for high in noTMS vs. lspcs
fprintf('Greater error for noTMS high vs lspcs high: %d\n',...
    sum(data_conds_noTMS_lspcs_high(:,1) - data_conds_noTMS_lspcs_high(:,2) > 0))

% how many subjects show greater error for low in noTMS vs. lspcs
fprintf('Greater error for noTMS low vs lspcs low: %d\n',...
    sum(data_conds_noTMS_lspcs_low(:,1) - data_conds_noTMS_lspcs_low(:,2) > 0))

%high no TMS vs lspcs; low no TMS vs. lspcs
[t,stats,orig,param] = permtest_t(data_noTMS_RVF_mat,data_lspcs_RVF_mat)



%get LVF descriptive stats
data_LVFW = data_aveW(~strcmp(data_aveW.condition,'l_ips2') & ...
    strcmp(data_aveW.hemi,'LVF'),:);

data_LVFW_ave = summary_stats(data_LVFW,{'mean_i_sacc_err_diff','mean_i_sacc_rt_diff'},{'condition','hemi'});


%% make efield figure (Fig 3)
efield_data = readtable(fullfile(data_dir,'efield_data_23-Feb-2024.csv'));

%wide version for stats
efield_dataW = unstack(efield_data,{'efield','voxels'},'hemi');

%average for plotting
efield_data.hemi01 = double(strcmp(efield_data.hemi,'rh'));
efield_data_ave = summary_stats(efield_data,'efield',{'roi','hemi01'})

f=figure('units','centimeters','Position',[20 20 10 12]);
hold on;
for ss=1:length(subj_valid)
    this_efield = efield_data(efield_data.subject==subj_valid(ss),:);
    plot(this_efield.hemi01,this_efield.efield, ...
        '-','Color',[0 0 0 0.5],'LineWidth',1);
end
errorbar(efield_data_ave.hemi01,efield_data_ave.mean_efield,efield_data_ave.sem_efield,'k')
xticks([0 1]);
xticklabels({'Left','Right'})
xlim([-0.5 1.5])
ylabel('Modeled electrical field (V/m)')

%save figure
exportgraphics(f,fullfile(fig_dir,sprintf('spcs_efield_%s.pdf', ...
    date())), ...
    'ContentType','vector')

% stats
[t,stats,orig,param] = permtest_t(efield_dataW.efield_lh,...
    efield_dataW.efield_rh)

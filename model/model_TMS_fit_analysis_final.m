% This script produces figures/analyses based on variable-precision model
% fits. Produces panels for Fig 3.
%
% Hallenbeck, Tardiff, et al., JNeurosci, 2025


clear all;close all
addpath('./model');


%% setup
%matlab defaults for appearance
set(groot,'DefaultFigureColor',[1 1 1])

set(groot,'DefaultAxesFontName','Arial')
set(groot,'DefaultLegendFontName','Arial');
set(groot,'DefaultTextFontName','Arial');
set(groot,'DefaultAxesFontSize',7); %15)

set(groot,'DefaultLineLineWidth',1.75)
set(groot,'DefaultErrorBarLineWidth',1.75)
set(groot,'DefaultErrorBarMarkerSize',12)
set(groot,'DefaultErrorBarMarkerFaceColor','auto')
set(groot,'DefaultErrorBarCapSize',0)
set(groot,'DefaultScatterLineWidth',1.75)
set(groot,'DefaultScatterMarkerFaceColor',[1 1 1]);
set(groot,'DefaultLegendBox','off')


%specify directories,conditions,subjects
fit_dir = './fits_nathan';
fig_dir = '../figs_raw';
data_dir = './data_061322';

valid_subj = [1 2 3 4 5 6 7 11 12 16 17 21 22 24]';

exptype = 'initial'; % saccade type: 'initial' or 'final'
condVec = {'noTMS','l_spcs'};
condVec_names = {'none','lSPCS'};
if iscell(condVec)
    condition = strjoin(condVec,'_');
else
    condition = condVec;
end
exppriorityVec = [2/3 1/3];
pconds = {'high','low'};
models = {'TMS_pJ'}; 

nLLVar = []; %'nLLVec_all';

e_max = 25;
de = 0.025;
erange = 0:de:e_max;
nerange = length(erange);

model_dirs = strcat(exptype,'_model_',models,'_',condition,'_subH');

%{
model_color = {[0,0.447,0.741]
                [0.85,0.325,0.098]
                [0.929,0.694,0.125]
                [0.494,0.184,0.556]
                [0.466,0.674,0.188]
                [0.301,0.745,0.933]
                [0.635,0.0780,0.184]};
%}

cond_color = lines(2);
if length(models)==1
    model_color = {[0 0 0]};
else
    model_color = lines(length(models));
    model_color = mat2cell(model_color,ones(1,length(models)));
end
            
%% load model fits and data 
params = struct();
fiteval_raw = [];
for m=1:length(models)
    [this_params,this_fiteval] = load_params(fullfile(fit_dir,model_dirs{m}),[],models{m},condVec,'all',nLLVar,valid_subj);

    %do any model-specific calculations here 
    if strcmp(models{m},'TMS_Jlow')
        this_params = compute_TMS_Jlow_params(this_params,condVec);
    end

    params.(models{m}) = this_params;
    fiteval_raw = [fiteval_raw; this_fiteval];
end

%get subj info
subjidVec = unique(fiteval_raw.subject);
nSubjs = length(subjidVec);

%make sure we have the correct subjects
assert(isequal(subjidVec,valid_subj),"Missing or invalid subjects detected!")


data_pred = [];
data_pred_ave = [];
for isubj = 1:nSubjs
    subjid = subjidVec(isubj);

    % load fitting data
    F=load(fullfile(data_dir,sprintf('data_subjid%d.mat',subjid)));
    this_data = [];
    this_data_ave = [];
    for c=1:length(condVec)
        this_cond = condVec{c};
        this_datac = cell2table([F.data.(exptype).(this_cond)',pconds'],'VariableNames',{'data','priority'});
        this_datac.condition = repmat(condVec(c),height(this_datac),1);

        for m=1:length(models)
            this_mod = models{m};
            if strcmp(this_mod,'jointsingle')
                this_paramflag = this_mod;
            else
                this_paramflag = [];
            end
            
            %average predictions across params
            this_params_raw = params.(this_mod)(params.(this_mod).subject==subjid,:);

            preds_raw = nan(nerange,length(condVec),height(this_params_raw));
            preds_ave_raw = nan(height(this_params_raw),2);
            preds_var_raw = nan(height(this_params_raw),2);
            for it=1:height(this_params_raw)
                this_params = get_param_values(this_params_raw(it,:),this_cond);
                [preds_raw(:,:,it),~,preds_ave_raw(it,:),preds_var_raw(it,:)] = ...
                    gen_preds(this_params,exppriorityVec,e_max,de,this_paramflag);
            end

            %average predictions for plotting/saving
            preds = mean(preds_raw,3);
            preds_ave = mean(preds_ave_raw,1);
            preds_var = mean(preds_var_raw,1);

           this_datac.(['pred_' this_mod]) = cell(2,1);
           for p=1:length(pconds)
                this_datac(strcmp(this_datac.priority,pconds{p}),['pred_' this_mod]) = ...
                    {preds(:,p)};
           end

           this_mdata_ave = this_datac(:,{'priority','condition'});
           this_mdata_ave.model = repmat({this_mod},height(this_mdata_ave),1);
           this_mdata_ave.pred_ave = preds_ave';
           this_mdata_ave.pred_var = preds_var';

           this_data_ave = [this_data_ave;this_mdata_ave];
        end

        
        this_data = [this_data;this_datac];
    end
    this_data.subject = repmat(subjid,height(this_data),1);
    this_data_ave.subject = repmat(subjid,height(this_data_ave),1);

    data_pred = [data_pred;this_data];
    data_pred_ave = [data_pred_ave;this_data_ave];

end

%add flags for priority (might not be using?)
data_pred.isHigh = double(strcmp(data_pred.priority,'high'));
data_pred_ave.isHigh = double(strcmp(data_pred_ave.priority,'high'));

%% plot model predictions against data


%summarise data for plotting
data_pred.data_ave = cellfun(@nanmean,data_pred.data);
data_pred.data_var = cellfun(@nanvar,data_pred.data);
data_ave_all = summary_stats(data_pred,{'data_ave','data_var'},{'condition','isHigh'});
data_ave_all.isHighC = categorical(data_ave_all.isHigh,[1,0],{'high','low'});
data_ave_all.conditionC = categorical(data_ave_all.condition,condVec);

% summarise model preds for plotting
data_pred_ave_all = summary_stats(data_pred_ave,{'pred_ave','pred_var'},{'condition','model','isHigh'});
data_pred_ave_all.isHighC = categorical(data_pred_ave_all.isHigh,[1,0],{'high','low'});
data_pred_ave_all.conditionC = categorical(data_pred_ave_all.condition,condVec);


%plot model means against data
f=figure('units','centimeters','Position',[20 20 12 12]);
hold on;
for c=1:length(condVec)
    %first plot effect in each condition
    this_data_ave = data_ave_all(data_ave_all.conditionC==condVec{c},:);
    errorbar(this_data_ave.isHighC,this_data_ave.mean_data_ave,...
        this_data_ave.sem_data_ave,'o','Color',cond_color(c,:));
    
    
    for v=1:length(models)
        this_pred_ave = data_pred_ave_all(data_pred_ave_all.conditionC==condVec{c} & ...
            strcmp(data_pred_ave_all.model,models{v}),:);
        plot(this_pred_ave.isHighC,this_pred_ave.mean_pred_ave,...
            '-','Color',cond_color(c,:));
            %'-','Color',model_color{v},'LineWidth',1);
    end
    
end
suplabel('priority','x',[0.0400 0.1100 0.9750 0.9200]); %,[0.0900 0.1100 0.4147 0.8950]); %,[0.0900 0.1100 0.8550 0.8950]);
ylabel('Error (\circ)')
ylim([1.4 2.4])
yticks([1.4:.5:2.4])
exportgraphics(f,fullfile(fig_dir,sprintf('TMS_pJ_fit_%s.pdf', date())), ...
    'ContentType','vector')


%% parameter analysis


%Jbar allocated to high and low
for c=1:length(condVec)
    for p=1:length(pconds)
        this_p = params.TMS_pJ.(['p_high_' condVec{c}]);
        if strcmp(pconds(p),'low')
            this_p = 1-this_p;
        end
        params.TMS_pJ.(['Jbar_' condVec{c} '_' pconds{p}]) = ...
            this_p.*params.TMS_pJ.(['Jbar_' condVec{c}]);
    end
end

%average params across iterations within subjects
params_TMS_pJ_ave_subj = summary_stats(params.TMS_pJ,...
    {'Jbar_noTMS','Jbar_l_spcs','tau','p_high_noTMS','p_high_l_spcs'...,
        'Jbar_noTMS_high','Jbar_noTMS_low','Jbar_l_spcs_high','Jbar_l_spcs_low'},...
        {'subject','model'});

%grand averages
params_TMS_pJ_ave = summary_stats(params_TMS_pJ_ave_subj,...
    strcat('mean_',{'Jbar_noTMS','Jbar_l_spcs','tau','p_high_noTMS','p_high_l_spcs'...,
        'Jbar_noTMS_high','Jbar_noTMS_low','Jbar_l_spcs_high','Jbar_l_spcs_low'}),...
        {'model'});
params_TMS_pJ_ave_med = varfun(@median,params.TMS_pJ,'InputVariables',...
    {'Jbar_noTMS','Jbar_l_spcs','tau','p_high_noTMS','p_high_l_spcs'...,
        'Jbar_noTMS_high','Jbar_noTMS_low','Jbar_l_spcs_high','Jbar_l_spcs_low'},...
        'GroupingVariables',{'model'});


%allocation p barplot-like figure for paper 
for c=1:length(condVec)
    params_TMS_pJ_ave.(['mean_mean_p_low_' condVec{c}]) = 1 - params_TMS_pJ_ave.(['mean_mean_p_high_' condVec{c}]);
end
p_aves = params_TMS_pJ_ave(:,startsWith(params_TMS_pJ_ave.Properties.VariableNames,'mean_mean_p'));
p_aves_mat = reshape(p_aves{:,:},2,2);

barpos = [0.4,0.25]; %explicitly setting this but should be this by default
barpos_offset = [-0.015,0.015];
bar_colors = viridis(2);
%subplot = @(m,n,p) subtightplot(m,n,p,.08,[0.11,0.05],[.08,0.025]);
f=figure('units','centimeters','Position',[20 20 16 12]);
%subplot(1,1,1)
hold on;

%plot markers
xline(0.5,':k','Alpha',1.0,'LineWidth',2)
xline(exppriorityVec(1),'--k','objective probability','Alpha',1.0,...
    'LabelVerticalAlignment','bottom','LabelOrientation','horizontal',...
    'LineWidth',1.5)

%plot bars
hb=barh(p_aves_mat,'stacked','XData',barpos,'FaceAlpha',0.4,'LineWidth',1.75,...
    'BarWidth',1);
hb(1).FaceColor = bar_colors(1,:);
hb(2).FaceColor = bar_colors(end,:);

%plot subject scatter
scatalpha = 0.75;
scatter(params_TMS_pJ_ave_subj{:,strcat('mean_p_high','_',condVec)}',...
    repmat(barpos+barpos_offset,height(params_TMS_pJ_ave_subj),1)',...
    80,'ko','MarkerEdgeAlpha',scatalpha,'MarkerFaceAlpha',scatalpha)
plot(params_TMS_pJ_ave_subj{:,strcat('mean_p_high','_',condVec)}',...
    repmat(barpos+barpos_offset,height(params_TMS_pJ_ave_subj),1)',...
    '-','Color',[0 0 0 scatalpha]);%,'LineWidth',1.5)

%plot errorbars
errorbar(params_TMS_pJ_ave{:,strcat('mean_mean_p_high_',condVec)},barpos,...
    [],[],...
    params_TMS_pJ_ave{:,strcat('sem_mean_p_high_',condVec)},...
    params_TMS_pJ_ave{:,strcat('sem_mean_p_high_',condVec)},...
    'k.','MarkerSize',1,'CapSize',0,'LineWidth',3,'Color',[0 0 0])

%ylim([0.3 0.35])
yticks(fliplr(barpos))
yticklabels(fliplr(condVec_names))
ylabel('TMS condition')
xticks([0,0.5,1])
xlabel('fit allocation')

exportgraphics(f,fullfile(fig_dir,sprintf('TMS_pJ_fit_allocation_%s.pdf', date())), ...
    'ContentType','vector')


%vertical verzion of Jbar(/tau) figure
%plot bars
barpos = [0.25,0.4]; %explicitly setting this but should be this by default

subplot = @(m,n,p) subtightplot(m,n,p,.08,[0.11,0.05],[.08,0.025]);
f=figure('units','centimeters','Position',[20 20 16 12]);
subplot(1,2,1);hold on;
hb = bar(params_TMS_pJ_ave{:,{'mean_mean_Jbar_noTMS','mean_mean_Jbar_l_spcs'}},'XData',barpos,...
    'LineWidth',1.75,'FaceColor','flat','FaceAlpha',0.4);
hb.CData = cond_color;

%plot subject scatter
plot(repmat(barpos,height(params_TMS_pJ_ave_subj),1)',...
    params_TMS_pJ_ave_subj{:,strcat('mean_Jbar_',condVec)}',...    
    '-','Color',[0 0 0 0.5]); %,'LineWidth',1.5)
%plot errorbars
errorbar(barpos,params_TMS_pJ_ave{:,strcat('mean_mean_Jbar_',condVec)},...
    params_TMS_pJ_ave{:,strcat('sem_mean_Jbar_',condVec)},...
    'k.','MarkerSize',1,'CapSize',0,'LineWidth',2)

yticks([0:3])
xticks(barpos);
xticklabels(condVec_names)
ylabel('fit capacity')

%now plot tau
subplot(1,2,2);hold on;
swarmchart(zeros(height(params_TMS_pJ_ave_subj),1),params_TMS_pJ_ave_subj.mean_tau,...
    'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5,'XJitterWidth',0.1);
errorbar(0,params_TMS_pJ_ave.mean_mean_tau,params_TMS_pJ_ave.sem_mean_tau,'-ok')
xticks([])
xlim([-0.1 .1])
ylabel('\tau')
yticks([0:.1:.4])
clear subplot

exportgraphics(f,fullfile(fig_dir,sprintf('TMS_pJ_Jbar_tau_%s.pdf', date())), ...
    'ContentType','vector')



%% plot schematics for relationship between Jbar and error rates/precision (empirical)
% distributions for paper
Jbar_cond_params = cell(length(condVec),length(pconds));
for c=1:length(condVec)
    for p=1:length(pconds)
        Jbar_cond_params{c,p} = ...
            sprintf('mean_Jbar_%s_%s',condVec{c},pconds{p});
    end
end
Jbar_min = min(params_TMS_pJ_ave_subj{:,Jbar_cond_params(:)},[],'all');
Jbar_max = max(params_TMS_pJ_ave_subj{:,Jbar_cond_params(:)},[],'all');

%choosing bounds based on empirical min/max Jbar we see in data
Jbar_min_plot = fix(Jbar_min*10)/10;
Jbar_max_plot = ceil(Jbar_max*10)/10;

Jbar_range = linspace(Jbar_min_plot,Jbar_max_plot,100);
exp_error = nan(length(Jbar_range),1);
for j=1:length(Jbar_range)
    %a bit of a hacky use of gen_preds, but it works.
    [~,~,exp_error(j),~] = ...
        gen_preds([Jbar_range(j),params_TMS_pJ_ave.mean_mean_tau,1],1,40,0.02);
end

J_range = linspace(0,3,100);
pcond_ls = {'-','--'};

subplot = @(m,n,p) subtightplot(m,n,p,.11,[0.11,0.05],[.08,0.025]);
f = figure('units','centimeters','Position',[20 20 16 12]);
subplot(2,1,1);hold on;
plot(Jbar_range,exp_error,'k');
xlabel('$\bar{J}$','Interpreter','latex')
ylabel('Error (\circ)')
ylim([0 5])

subplot(2,1,2);hold on;
xlabel('J');
ylabel('Probability Density')

for c=1:length(condVec)
    for p=1:length(pconds)
        this_J_param_name = sprintf('mean_%s',Jbar_cond_params{c,p});
        this_J_param = params_TMS_pJ_ave.(this_J_param_name);

        subplot(2,1,1);hold on;
        xline(this_J_param,pcond_ls{p}, ...
            'Color',cond_color(c,:),'LineWidth',1.5,'Alpha',1);

        subplot(2,1,2);hold on;
        this_pJ = gampdf(J_range,...
            this_J_param/params_TMS_pJ_ave.mean_mean_tau,...
            params_TMS_pJ_ave.mean_mean_tau); % probability of that J value
        plot(J_range,this_pJ,pcond_ls{p},'Color',cond_color(c,:));
        
    end
end
exportgraphics(f,fullfile(fig_dir,sprintf('VP_precision_%s.pdf', ...
    date())), ...
    'ContentType','vector')
clear subplot



%% hypothesis schematics

tau_schem = params_TMS_pJ_ave.mean_mean_tau;
anoise = 0;
Jbar_storage = [1.6 1.15];
%keep mean error based on jbar fixed as alpha changes
%Jbar_storage = Jbar_noise_adj(Jbar_storage,anoise); %Jbar_storage./(1-anoise.*sqrt(Jbar_storage)).^2;
p_priority_noTMS = exppriorityVec;
p_high_TMS = [p_priority_noTMS(1) 0.5];


% H1: disrupt storage
f=figure('units','centimeters','Position',[20 20 16 16]);
%subplot(1,2,1);hold on;
hold on;
for c=1:length(condVec)
    [~,~,this_pred] = ...
        gen_preds([Jbar_storage(c),tau_schem,anoise,p_priority_noTMS(1)],p_priority_noTMS,40,0.02,'anoise');
        %gen_preds([Jbar_storage(c),tau_schem,p_priority_noTMS(1)],p_priority_noTMS,40,0.02);

    plot(categorical(pconds),this_pred,...
        '-');%,'Color',cond_color(c,:));
    
end
%ylabel('Error (\circ)')
%title('disrupt storage')
%yticks([])
%ylim([0 3])
%suplabel('priority','x',[0.0400 0.1100 0.9750 0.9200]); %,[0.0900 0.1100 0.4147 0.8950]); %,[0.0900 0.1100 0.8550 0.8950]);

% H2: disrupt prioritization
%let's actually do this over a range of Jbars
Jbar_priority = 1.6;
%f=figure('units','centimeters','Position',[20 20 16 16]);
%subplot(1,2,2)
%hold on;
for c=2 %1:length(condVec)
    [~,~,this_pred] = ...
        gen_preds([Jbar_priority,tau_schem,anoise,p_high_TMS(c)],p_priority_noTMS,40,0.02,'anoise');

    plot(categorical(pconds),this_pred,...
        '-');%,'Color',cond_color(c,:));
    
end
ylim([0 3])
%ylim([1.1 2.6])
%suplabel('priority','x',[0.0400 0.1100 0.9750 0.9200]); %,[0.0900 0.1100 0.4147 0.8950]); %,[0.0900 0.1100 0.8550 0.8950]);
ylabel('Error (\circ)')
yticks([])
%legend(condVec,'Location','northwest')
%title('disrupt prioritization')
%title('disrupt priority')
suplabel('priority','x',[0.0400 0.1100 0.9750 0.9200]); %,[0.0900 0.1100 0.4147 0.8950]); %,[0.0900 0.1100 0.8550 0.8950]);

ll=legend({'none','disrupt storage','disrupt priority'},'Location','northwest');
title(ll,'TMS effect')
exportgraphics(f,fullfile(fig_dir,sprintf('VP_model_preds_%s.pdf', ...
    date())), ...
    'ContentType','vector')


% now plot the error vs. Jbar predictions for each model

subplot = @(m,n,p) subtightplot(m,n,p,.11,[0.11,0.05],[.08,0.025]);
f = figure('units','centimeters','Position',[20 20 16 12]);

%now plot H1 Jbar/predictions
subplot(2,1,1);hold on;
plot(Jbar_range,exp_error,'k');
xlabel('$\bar{J}$','Interpreter','latex')
ylabel('Error (\circ)')
ylim([0 5])
for c=1:length(condVec)
    for p=1:length(p_priority_noTMS)
        this_bar = Jbar_storage(c)*p_priority_noTMS(p);
        xline(this_bar,pcond_ls{p},'Color',cond_color(c,:));
    end
end
%set(gca,'Xscale','log')

%now plot H2 Jbar/predictions
subplot(2,1,2);hold on;
plot(Jbar_range,exp_error,'k');
xlabel('$\bar{J}$','Interpreter','latex')
ylabel('Error (\circ)')
ylim([0 5])

for c=1:length(condVec)
    this_pvec = [p_high_TMS(c) 1-p_high_TMS(c)];
    for p=1:length(this_pvec)
        this_bar = Jbar_priority*this_pvec(p);
        xline(this_bar,pcond_ls{p},'Color',cond_color(c,:));
    end
end
%set(gca,'Xscale','log')

clear subplot

exportgraphics(f,fullfile(fig_dir,sprintf('VP_error_J_preds_%s.pdf', ...
    date())), ...
    'ContentType','vector')


%% statistical testing

%do we replicate Yoo et al. 2018 finding of priority effect?
%yes, but this was guarnteed by the subject inclusion procedure, so
%probably don't need to report, unless as a validation of the model
[t stats orig param] = permtest_t(params_TMS_pJ_ave_subj.mean_p_high_noTMS,[],'m',0.5)
stats.ci + 0.5
fprintf('p_high ~= 0.5 statistically in noTMS group \n')


%do we replicate Yoo et al. 2018 finding of under-allocating
%to high probability without TMS?
[t stats orig param] = permtest_t(params_TMS_pJ_ave_subj.mean_p_high_noTMS,[],'m',exppriorityVec(1))
stats.ci + exppriorityVec(1)

fprintf('p_high < objective p_high in %d subjects in no TMS and %d in l_spcs out of %d \n',...
    sum(params_TMS_pJ_ave_subj.mean_p_high_noTMS<exppriorityVec(1)),...
    sum(params_TMS_pJ_ave_subj.mean_p_high_l_spcs<exppriorityVec(1)),length(subjidVec))


%is there a difference in p_high between conds? (yes!)
[t stats orig param] = permtest_t(params_TMS_pJ_ave_subj.mean_p_high_noTMS,params_TMS_pJ_ave_subj.mean_p_high_l_spcs)
stats.ci

%number w p greater for no TMS?
fprintf('p_high greater in no TMS in %d out of %d subjects\n',...
    sum(params_TMS_pJ_ave_subj.mean_p_high_noTMS>params_TMS_pJ_ave_subj.mean_p_high_l_spcs),length(subjidVec))


%is there a difference in Jbar between conds? (yes, weakly)
%perm test version
[t stats orig param] = permtest_t(params_TMS_pJ_ave_subj.mean_Jbar_noTMS,params_TMS_pJ_ave_subj.mean_Jbar_l_spcs)
stats.ci

%number w/ Jbar lower for no TMS?
fprintf('Jbar lower in no TMS in %d out of %d subjects\n',...
    sum(params_TMS_pJ_ave_subj.mean_Jbar_noTMS<params_TMS_pJ_ave_subj.mean_Jbar_l_spcs),length(subjidVec))
%mean Jbar by condition
fprintf('Jbar noTMS: %.2f; Jbar l_spcs: %.2f; diff: %.2f\n',...
    params_TMS_pJ_ave.mean_mean_Jbar_noTMS, params_TMS_pJ_ave.mean_mean_Jbar_l_spcs,...
    params_TMS_pJ_ave.mean_mean_Jbar_noTMS-params_TMS_pJ_ave.mean_mean_Jbar_l_spcs)

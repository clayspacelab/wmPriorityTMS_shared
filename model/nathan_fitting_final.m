% Code for fitting the variable-precision priority model to data
%
% Hallenbeck, Tardiff, et al., JNeurosci, 2025


clear all
addpath('./model','helper_functions','-end');

data_dir = 'data_061322'; 
fit_dir = 'fits_nathan';
overwrite = false;

exptype = 'initial'; % saccade type: 'initial' or 'final'
model = 'TMS_pJ'; 
condVec = {'noTMS','l_spcs'};
subsampHigh = true; %subsample the high priority trials
nConds = length(condVec);

% subject list
subjidVec = [1 2 3 4 5 6 7 11 12 16 17 21 22 24];
nSubjs = length(subjidVec);

runmax = 20; % how many fits to run per participant
runlist = 1:runmax;
exppriorityVec = [2/3 1/3]; % experimental priority vector
fixparams = []; 
nfixparams = size(fixparams,2);

if subsampHigh
    warning('nsamps for BIC calc will not be correct for subsampled data!')
end
for isubj = 1:nSubjs
    subjid = subjidVec(isubj);

    % load fitting data
    load(fullfile(data_dir,sprintf('data_subjid%d.mat',subjid)))
    if iscell(condVec)
        data = data.(exptype);
    else
        data = data.(exptype).(condVec);
    end

    %compute total sample size for later calc of BIC
    if iscell(condVec)
        nsamps = 0;
        for c = 1:length(condVec)
            nsamps = nsamps + sum(cellfun(@(x) sum(~isnan(x)),data.(condVec{c})));
        end
    else
        nsamps = sum(cellfun(@(x) sum(~isnan(x)),data));
    end
    
    % file saving name
    if iscell(condVec)
        condition = strjoin(condVec,'_');
    else
        condition = condVec;
    end
    filepath = fullfile(fit_dir,sprintf('%s_model_%s_%s',exptype,model,condition));
    if subsampHigh
        filepath = [filepath '_subH'];
    end
    filename = sprintf('%s/fits_subj%d.mat',filepath,subjid);
    if ~exist(filepath,'dir'); mkdir(filepath); end
    if ~overwrite && exist(filename,'file')
        error('File already exists!')
    end
    
    ML_parameters = []; nLLVec = [];
    runlist_completed = [];
    if subsampHigh
        nLLVec_all = [];
    end

    %try load(filename); catch; ML_parameters = []; nLLVec = []; end
    %try runlist_completed*2; catch; runlist_completed = []; end % seeing if runlist_completed exists yet

    for irun = 1:length(runlist)
        runlistt = runlist(irun);
        rng(runlistt)
        %try
           irun
            subjid
            %NOTE: fit_parameters has its own inner loop, but for whatever
            %reason this code is running it in an outer loop and feeding in
            %iterations (with iteration numbers) one at a time. This probably slows things down a
            %bit but not worth fixing atm.
            if subsampHigh
                fit_data = do_subsampH(data,subjid,irun,condVec);
            else
                fit_data = data;
            end

            [bfp,fval,rlc] = fit_parameters(model,fit_data,exppriorityVec,runlistt,runmax,fixparams,condVec);
            
            ML_parameters = [ML_parameters; bfp];
            nLLVec = [nLLVec fval];
            runlist_completed = [runlist_completed rlc];

            if subsampHigh
                fval_all = calc_nLL_postSub(model,data,bfp,exppriorityVec,condVec);
                nLLVec_all = [nLLVec_all fval_all];
                save(filename,'ML_parameters','nLLVec','nLLVec_all','runlist_completed','nsamps','condVec','nfixparams')
            else
                save(filename,'ML_parameters','nLLVec','runlist_completed','nsamps','condVec','nfixparams')
            end
        %end
    end
end
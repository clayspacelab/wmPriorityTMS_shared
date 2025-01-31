function [ML_parameters, nLLVec, runlist_completed] = fit_tms_stronghyp(data,exppriorityVec,runlist,runmax,fixparams,condhypVec)
%FIT_PARAMETERS estimates the parameters given model and data
% 
%   FIT_PARAMETERS(MODEL, DATA, EXPPRIORITYVEC) estimates maximum-likelihood 
%     parameters for MODEL, DATA, and EXPPRIORITYVEC.
% 
%   FIT_PARAMETERS(MODEL, DATA, EXPPRIORITYVEC, RUNLIST, RUNMAX) estimates 
%     maximum-likelihood parameters for MODEL, DATA, and EXPPRIORITYVEC for
%     a vector of indices RUNLIST contained in RUNMAX.
%
%   FIT_PARAMETERS(MODEL, DATA, EXPPRIORITYVEC, RUNLIST, RUNMAX, FIXPARAMS)
%     uses FIXPARAMS to fix select parameters to a specified value, taking
%     those parameters out of the optimization step.
% 
%   ============ INPUT VARIABLES ============
% 
%   MODEL: string. 'max_points', 'flexible', 'proportional', or 'min_error'
% 
%   DATA: cell of length nPriorities. the ith element of DATA should be
%     all data corresponding to EXPPRIORITYVEC(i) condition. the first
%     column should contain the magnitude of errors and the second column,
%     if available, should contain the corresponding circle wager radius
%     size. 
%
%   EXPPRIORITYVEC: vector of priority values (in decreasing order of
%     priority) e.g., [0.6 0.3 0.1]
% 
%   RUNLIST: (optional). scalar or vector containing indices 1 to RUNMAX. 
%     this variable is useful when parallelizing the MLE step.
% 
%   RUNMAX: (optional). total number of optimizations you plan to start for 
%     a given TESTMODEL and SUBJNUM. spaces out the starting values using
%     a latin hypercube. this uses a seed, necessary for replicability.
%     this variable is useful when parallelizing the MLE step.
% 
%   FIXPARAMS: (optional). 2 x (number of fixed parameters) matrix. fixed 
%     parameters, such that the first row corresponds to the index and 
%     second row corresponds to the value of the fixed parameter. 
%   
%   ========== OUTPUT VARIABLES ==========
%   
%   ML_PARAMETERS: length(runlist) x nParameters matrix. each row of this
%     matrix corresponds to the MLE of one optimization. default
%     length(runlist) == 1. 
% 
%   NLL_VEC: vector of length length(runlist), containing the negative
%     log-likelihoods of each of the length(runlist) runs. default
%     length(runlist) == 1. 
%
%   RUNLIST_COMPLETED: vector of indices, RUNLIST, completed. useful when
%     parallelizing the MLE step.

% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% -----------------------

if nargin < 4; runlist = 1; end
if nargin < 5; runmax = 50; end
if nargin < 6; fixparams = []; end

[logflag, lb, ub, plb, pub] = loadconstraints('stronghyp',exppriorityVec);

% VERY POORLY CODED: gtting rid of parietal condition
if sum(condhypVec == 0)
    switch (condhypVec(3))
        case 2
            fixparams = [4; 0];
        case 3
            fixparams = [5; 0];
    end
end

if ~(isempty(fixparams))
    % vector of non-fixed parameter indices
    nParams= length(logflag);
    freeparamsidx = 1:nParams;
    freeparamsidx(fixparams(1,:)) = [];
    
    logflag(fixparams(1,:)) = [];
    lb(fixparams(1,:)) = [];
    ub(fixparams(1,:)) = [];
    plb(fixparams(1,:)) = [];
    pub(fixparams(1,:)) = [];
end

% create list of all x0s
nParams = length(logflag);
x0_list = [];
constantt = 0;
while size(x0_list,1) < runmax
    rng(0);
    x0_list = lhs(runmax + constantt,nParams,plb,pub,[],1e3);
%     isviolated = check_nonbcon(model,x0_list);
%     x0_list(logical(isviolated),:) =[];
%     constantt = constantt + sum(isviolated);
end
x0_list(:,logflag) = exp(x0_list(:,logflag));

% optimize for starting values in RUNLIST
ML_parameters = [];
nLLVec = [];
runlist_completed = [];
for irun = 1:length(runlist)
    runlist(irun)
   
    rng(runlist(irun));
    
    x0 = x0_list(runlist(irun),:);
    fun = @(x) calc_nLL_stronghyp(x,data,exppriorityVec,fixparams,condhypVec);
    
    x = bads(fun,x0,lb,ub,plb,pub);
    fval = fun(x);
    
    x(logflag) = exp(x(logflag));
    if isempty(fixparams)
        bfp = x;
    else
        bfp = nan(1,nParams+size(fixparams,2));
        bfp(freeparamsidx) = x;
        bfp(fixparams(1,:)) = fixparams(2,:);
    end

    ML_parameters = [ML_parameters; bfp];
    nLLVec = [nLLVec fval];
    runlist_completed = [runlist_completed runlist(irun)];

end
end


function [logflag, lb, ub, plb, pub, nonbcon] = loadconstraints(model,exppriorityVec,is_wagerdata)
%LOADCONSTRAINTS loads optimization constraints for each model
% 
% LOADCONSTRAINTS(MODEL,EXPPRIORITYVEC) loads optimization constraint for a
% particular model
%
% LOADCONSTRAINTS(MODEL,EXPPRIORITYVEC,IS_WAGERDATA) loads optimization
% constraint for wager data, if IS_WAGERDATA = 1
% 
% ===== INPUT VARIABLES =====
% 
% MODEL: 'max_points', 'flexible', 'proportional', or 'min_error'
% 'flex_jbars'
% EXPPRIORITYVEC: row vector of experimental priority.
%   sum(exppriorityVec) = 1
%
% IS_WAGERDATA: 0: just error data. 1: error and wager data

if nargin < 3; is_wagerdata = 0; end

% change values if fitting all conditions jointly
if strcmp(model(1:3),'all')
    model = model(5:end); 
    allflag = 1;
else
    allflag = 0;
end

nPriorities = length(exppriorityVec); 

% lower and upper bounds, logflags (this is all that will be set if model
% is 'vp' (i.e. basic variable-precision model)
lb = [0.01 0.01]; % Jbar_total, tau
ub = [50 10];
plb = [0.5 0.01];
pub = [10 1];
logflag = [1 1];

if (allflag) % if fitting all TMS conditions jointly
    lb = lb([1 1 2]); % Jbar_total, Jbar_total, tau
    ub = ub([1 1 2]);
    plb = plb([1 1 2]);
    pub = pub([1 1 2]);
    logflag = logflag([1 1 2]);
end

if (is_wagerdata) % alpha beta
    lb = [lb 1e-5 1e-5];
    ub = [ub 5 5];
    plb = [plb 0.7 0.5];
    pub = [pub 1.3 1.5];
    logflag = [logflag 0 0];
end

if endsWith(model,'anoise')
    lb = [lb 0];
    ub = [ub 1];
    plb = [plb 0.08];
    pub = [pub 0.6];
    logflag = [logflag 0];
end

switch model
    case {'flexible','stronghyp','base','TMS_pJ',...
            'flexible_anoise'} % define p_high p_med for Flexible model
        lb = [lb 1e-10.*ones(1,nPriorities-1)];
        ub = [ub ones(1,nPriorities-1)-eps];
        plb = [plb max([1e-10.*ones(1,nPriorities-1); exppriorityVec(1:end-1).*0.5])];
        pub = [pub min([ones(1,nPriorities-1)-eps; exppriorityVec(1:end-1).*1.5])];
        logflag = [logflag zeros(1,nPriorities-1)];
        
        if (allflag) % if fitting all TMS conditions jointly
            lb = [lb lb(end).*ones(1,nPriorities-1)]; 
            ub = [ub ub(end).*ones(1,nPriorities-1)];
            plb = [plb plb(end).*ones(1,nPriorities-1)];
            pub = [pub pub(end).*ones(1,nPriorities-1)];
            logflag = [logflag logflag(end).*ones(1,nPriorities-1)];
        end
        
    case 'min_error' % gamma for Minimizing Error model
        lb = [lb 1e-10];
        ub = [ub 5];
        plb = [plb 1e-3];
        pub = [pub 1];
        logflag = [logflag 1];
        nonbcon = @model4nonbcon; % violates if Jbar/tau - gamma/2 <= 0 
    otherwise
        nonbcon = [];
end

switch model
    case 'stronghyp' % adds another jbar and p's to fit
        lb = [lb 1e-5 1e-10.*ones(1,nPriorities-1)];
        ub = [ub 50 ones(1,nPriorities-1)-eps];
        plb = [plb 0.5 max([1e-10.*ones(1,nPriorities-1); exppriorityVec(1:end-1).*0.5])];
        pub = [pub 10 min([ones(1,nPriorities-1)-eps; exppriorityVec(1:end-1).*1.5])];
        logflag = [logflag 1 zeros(1,nPriorities-1)];
end

switch model
    case 'flex_jbars'
        lb = [1e-5 1e-5 1e-3]; % Jbar_high, Jbar_low, tau
        ub = [50 50 10];
        plb = [0.5 0.5 0.01];
        pub = [10 10 1];
        logflag = [1 1 1];
    case 'both_flex_jbars'
        lb = [1e-5 1e-5 1e-5 1e-5 1e-3]; % Jbar_high, Jbar_low, tau
        ub = [50 50 50 50 10];
        plb = [0.5 0.5 0.5 0.5 0.01];
        pub = [10 10 10 10 1];
        logflag = [1 1 1 1 1];
end

%Nathan's models (not being cute and just having as many ps as there are in
%this study...)
%unfortunately I can't just define my own b/c the inner loop uses flexible
%model...THIS ASSUMES ONLY ONE p (i.e. two priorities) per TMS condition!
switch model
    case {'TMS_pJ'} %Jbar_noTMS Jbar_lspcs tau p_noTMS p_lspcs
        lb = [lb(1) lb lb(end)];
        ub = [ub(1) ub ub(end)];
        plb = [plb(1) plb plb(end)];
        pub = [pub(1) pub pub(end)];
        logflag = [logflag(1) logflag logflag(end)];
    case 'jointsingle' %Jbar/tau for each condition, order s/b (noTMS_high,noTMS_low,lspcs_high,lspcs_low)
        lb = repmat(lb,1,4);
        ub = repmat(ub,1,4);
        plb = repmat(plb,1,4);
        pub = repmat(pub,1,4);
        logflag = repmat(logflag,1,4);
end

logflag = logical(logflag); 
lb(logflag) = log(lb(logflag)); 
ub(logflag) = log(ub(logflag)); 
plb(logflag) = log(plb(logflag)); 
pub(logflag) = log(pub(logflag)); 


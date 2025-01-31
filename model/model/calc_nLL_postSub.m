function [nLL] = calc_nLL_postSub(model,data,params,exppriorityVec,condVec)
%calc_nLL_postSub computes nLL over all data for a given model. Useful if
%model was originally fit on subsampled data.
% 
%
% Nathan Tardiff. Adapted from code by Aspen Yoo.

if ~exist('condVec','var')
    condVec = {};
end

if ~(strcmp(model(1:3),'all')) && ~(strcmp(model(1:4),'both')) && ...
        ~strcmp(model,'TMS_p') && ~strcmp(model,'TMS_J') && ~strcmp(model,'TMS_Jlow') ...
        && ~strcmp(model,'base') && ~strcmp(model,'jointsingle') && ~strcmp(model,'TMS_pJ') ...
        && ~strcmp(model,'TMS_pJ_anoise')
    expnumber = size(data{1},2);
    assert(length(data) == length(exppriorityVec));
    if strcmp(model,'max_points')
        assert(expnumber == 2, 'maximizing points model only works with wager data')
    end
else
    expnumber = 1;
end

fixparams = [];


%take log of params that were fit on log scale
logflag = loadconstraints(model,exppriorityVec,expnumber-1);
x = params;
x(logflag) = log(x(logflag));
    
if strcmp(model(1:3),'all') 
    fun = @(x) calc_nLL_allcond(model,x,data,exppriorityVec,fixparams);
elseif strcmp(model(1:4),'both')
    fun = @(x) calc_nLL_noTMS_l_spcs(model,x,data,exppriorityVec,fixparams);
elseif strcmp(model,'base')
    fun = @(x) calc_nLL_base(x,data,exppriorityVec,fixparams,condVec);
elseif strcmp(model,'jointsingle')
    fun = @(x) calc_nLL_jointsingle(x,data,exppriorityVec,fixparams,condVec);
elseif strcmp(model,'TMS_p')
    fun = @(x) calc_nLL_TMS_p(x,data,exppriorityVec,fixparams,condVec);
elseif strcmp(model,'TMS_J')
    fun = @(x) calc_nLL_TMS_J(x,data,exppriorityVec,fixparams,condVec);
elseif strcmp(model,'TMS_pJ')
    fun = @(x) calc_nLL_TMS_pJ(x,data,exppriorityVec,fixparams,condVec);
elseif strcmp(model,'TMS_pJ_anoise')
    fun = @(x) calc_nLL_TMS_pJ_anoise(x,data,exppriorityVec,fixparams,condVec);
elseif strcmp(model,'TMS_Jlow')
    fun = @(x) calc_nLL_TMS_Jlow(x,data,exppriorityVec,fixparams,condVec,logflag);
else
    fun = @(x) calc_nLL(model,x,data,exppriorityVec,fixparams);
end

%compute nLL for given model
nLL = fun(x);

end


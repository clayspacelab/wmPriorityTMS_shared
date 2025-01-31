function [params,fiteval] = load_params(fit_dir,param_names,model_name,condVec,type,nLLVar,subjlist)%,save_fits)
%loads and saves best fit parameter values and fit stats (e.g. NLL, AIC, BIC) in
%tables

%not implemented yet!
%{
if ~exist('save_fits','var')
    save_fits = true;
end
%}

%add the ability to return only for best fit iteration, or return
%everything
if ~exist('type','var') || isempty(type)
    type = 'best';
end

%allows the ability to specify the name of the nLL variable, which is used
%when computing fit stats using all data for fits fit on subsampled data.
if ~exist('nLLVar','var') || isempty(nLLVar)
    nLLVar = 'nLLVec';
end

fit_files = struct2table(dir(fullfile(fit_dir,"fits_subj*.mat")));

subjnum = regexp(fit_files.name,'fits_subj(\d+)','tokens');
fit_files.subjnum = cellfun(@(x) str2double(x{1}),subjnum);

if exist('subjlist','var') && ~isempty(subjlist)
    fit_files(~ismember(fit_files.subjnum,subjlist),:) = [];
end

if isempty(param_names)
    param_names = get_param_names(char(model_name),condVec);
end

if ischar(model_name)
    model_name = cellstr(model_name);
end

params = [];
fiteval = [];
for i=1:height(fit_files)
    F = load(fullfile(fit_files.folder{i},fit_files.name{i}));
    %find best fit
    if strcmp(type,'best')
        [~,mIdx] = min(F.(nLLVar));
        this_params = array2table(F.ML_parameters(mIdx,:),'VariableNames',param_names);
    elseif strcmp(type,'all')
        this_params = array2table(F.ML_parameters,'VariableNames',param_names);
    end
    this_params.subject = repmat(fit_files.subjnum(i),height(this_params),1);
    this_params.model = repmat(model_name,height(this_params),1);
    params = [params;this_params];

    if strcmp(type,'best')
        this_fiteval = table(F.(nLLVar)(mIdx),fit_files.subjnum(i),model_name,'VariableNames',{'nlle','subject','model'});
    elseif strcmp(type,'all')
        this_nlle = F.(nLLVar)';
        this_fiteval = table(this_nlle,repmat(fit_files.subjnum(i),length(this_nlle),1),...
            repmat(model_name,length(this_nlle),1),'VariableNames',{'nlle','subject','model'});
    end
    if isfield(F,'nsamps')
        this_fiteval.nsamps = repmat(F.nsamps,height(this_fiteval),1);
    end
    if isfield(F,'nfixparams')
        this_fiteval.nfixparams = repmat(F.nfixparams,height(this_fiteval),1);
        nfixflag = 0;
    else
        this_fiteval.nfixparams = zeros(height(this_fiteval),1);
        nfixflag = 1;
    end

    fiteval = [fiteval; this_fiteval];
end

%sort tables by subject
params = sortrows(params,'subject');
fiteval = sortrows(fiteval,'subject');

%extract number of fixed params if present
if nfixflag
    warning('No fixed parameter info detected. Defaulting to 0!')
end
assert(length(unique(fiteval.nfixparams))==1,'Cannot process files w/ variation in fixed params!')
nfixparams = fiteval.nfixparams(1);

%add AIC to fiteval (can't do BIC w/o info on number of samples, which I
%should start storing...)
k = length(param_names) - nfixparams;
fiteval.nparams = repmat(k,height(fiteval),1);
fiteval.AIC = 2.*k + 2.*fiteval.nlle;  %already negative loglik, so add not subtract
%add BIC if samples info available
if ismember('nsamps',fiteval.Properties.VariableNames)
    fiteval.BIC = log(fiteval.nsamps).*k + 2.*fiteval.nlle;
else
    warning('Number of samples not available for computing BIC!')
end


end
function [pred_data,erange,pred_data_ave,pred_data_var] = gen_preds(params,expPriorityVec,max_error,de,param2)

%we are going to devise the analytical error distributions for each
%priority level, and if requested by inclusion of nTrials, will simulate
%nTrials draws from the distribution (drawn proprortionaly to expPriorityVec).

%implemented for base model with one Jbar/tau, and some number of priority
%levels. Implement model-specific functions to handle more complex cases
%using this as the base function.

if ~exist('max_error','var') || isempty(max_error)
    max_error = 15; %maximum error to account for in distribution
end
if ~exist('de','var') || isempty(de)
    de = 0.02;
end
if ~exist('param2','var') || isempty(param2)
    param2 = 'scale';
end

nPriorities = length(expPriorityVec);
%nerange = 500;

%leaving this as array for backward compatibility with old fits
if ~strcmp(param2,'jointsingle')
    Jbar_total = params(1);
    tau = params(2);
    pp = params(end-(nPriorities-2):end);
    pVec = [pp 1-sum(pp)];
end
if strcmp(param2,'anoise')
    anoise = params(3);
else
    anoise = 0;
end

%erange = linspace(0,max_error,nerange); %error range for which to calculate probabilities.
erange = 0:de:max_error;
nerange = length(erange);

pred_data = nan(nerange,nPriorities);
for ipriority = 1:nPriorities
    if strcmp(param2,'jointsingle')
        %we have passed in individual Jbars and taus for each priority
        %rather than a priority with a Jtotal and a fixed tau.
        Jbar = params(1+(ipriority-1)*2);
        tau = params(2+(ipriority-1)*2);
    else
        Jbar = pVec(ipriority)*Jbar_total;
    end
    %{
    if ismember(param2,{'scale','jointsingle','anoise'})
        %default: tau is scale param, Jbar is mean
        [JVec] = loadvar('JVec',{Jbar,tau}); % values of J
        Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
    elseif strcmp(param2,'shape') 
        %in "shape" case, tau is actually shape not scale
        [JVec] = loadvar('JVec',{Jbar,Jbar/tau}); % values of J
        Jpdf = gampdf(JVec,tau,Jbar/tau); % probability of that J value
    end
    %}
    if strcmp(param2,'shape') 
        %in "shape" case, tau is actually shape not scale
        [JVec] = loadvar('JVec',{Jbar,Jbar/tau}); % values of J
        Jpdf = gampdf(JVec,tau,Jbar/tau); % probability of that J value
    else
        %default: tau is scale param, Jbar is mean
        [JVec] = loadvar('JVec',{Jbar,tau}); % values of J
        Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
    end
    Jpdf_int = trapz(JVec,Jpdf);
    if abs(Jpdf_int-1) > 1e-2
        warning('pdf does not integrate to 1: %f (J=%f,tau=%f,p#=%d)!',Jpdf_int,Jbar,tau,ipriority)
    end
    %assert(trapz(JVec,Jpdf)-1 < 1e-2,'pdf does not integrate to 1: %f (J=%f,tau=%f)!',trapz(JVec,Jpdf),Jbar,tau)
    Jpdf = Jpdf./Jpdf_int; % normalize slop

    nJs = length(JVec);
    p_Shat = raylpdf(repmat(erange,nJs,1),repmat(sqrt(1./JVec)'+anoise,1,length(erange)));
    p_Shat(p_Shat == 0) = 1e-10; % set to arbitrarily small value if zero

    pTrials = trapz(JVec',p_Shat.*Jpdf');
    pTrials_int = trapz(erange,pTrials);
    if abs(pTrials_int-1) > 1e-2
        warning('marginal error distribution does not integrate to 1: %f (J=%f,tau=%f,p#=%d)',pTrials_int,Jbar,tau,ipriority);
    end

    pTrials = pTrials./pTrials_int;    %fix slop
    
    %store marginal distribution by priority
    pred_data(:,ipriority) = pTrials;
end

%compute mean of distribution
pred_data_ave = trapz(erange,pred_data.*erange');

%compute variance of distribution
pred_data_var = trapz(erange,pred_data.*(erange'-pred_data_ave).^2);



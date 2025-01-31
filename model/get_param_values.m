function [param_vals] = get_param_values(params,cond,model)
%given a known model, return the VP model parameters associated with a
%given TMS condition (i.e. Jbar,tau,p_high)
%for now we assume params is a table or struct...could relax this at some
%point if needed

if ismember('model',params.Properties.VariableNames)
    model = char(params.model);
elseif ~exists('model','var')
    error('Must supply model name if not passed in params')
end

%condition-specific parameter names
p_high_cond = ['p_high_' cond];
Jbar_cond = ['Jbar_' cond];

switch model
    case {'flexible','base'}
        %does not take into account TMS condition
        param_names = {'Jbar','tau','p_high'};
    case {'TMS_pJ'}
        param_names = {Jbar_cond,'tau',p_high_cond};
    case 'jointsingle'
        param_names_base = {Jbar_cond ['tau_' cond]};
        param_names = [strcat(param_names_base,'_high'),strcat(param_names_base,'_low')];
end

param_vals = params{:,param_names};

end
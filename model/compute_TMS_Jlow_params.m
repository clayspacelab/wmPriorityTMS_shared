function [params] = compute_TMS_Jlow_params(raw_params,condVec)
%compute TMS_Jlow model derived parameters (Jbar, Jlow)
if ischar(condVec)
    condVec = {condVec};
end

params = raw_params;
for c=1:length(condVec)
    this_p = params.(['p_high_' condVec{c}]);
    params.(['Jlow_' condVec{c}]) = ...
        (params.Jhigh.*(1-this_p))./this_p;
    params.(['Jbar_' condVec{c}]) = ...
        params.Jhigh + params.(['Jlow_' condVec{c}]);
end

end
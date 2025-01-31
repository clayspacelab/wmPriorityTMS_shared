function [param_names] = get_param_names(model,condVec)

switch model
    case {'flexible','base'}
        param_names = {'Jbar','tau','p_high'};
    case 'TMS_pJ'
        param_names = {'Jbar','Jbar','tau','p_high','p_high'};
        %for some reason strcat won't operate when param names are
        %different sizes...??
        param_names(1:2) = strcat(param_names(1:2),'_',condVec);
        param_names(4:5) = strcat(param_names(4:5),'_',condVec);
    case 'jointsingle'
        param_names_base = {'Jbar','tau'};
        param_names = {};
        for c=condVec
            for p={'high','low'}
                param_names = [param_names strcat(param_names_base,...
                    '_',c,'_',p)];
            end
        end
end

end
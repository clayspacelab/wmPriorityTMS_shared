function nLL = calc_nLL_TMS_pJ(Theta,data,exppriorityVec,fixparams,condVec)
%effect of TMS is to change priority, not Jtotal

nLL = 0;

%NOTE THIS ASSUMES TWO CONDITIONS. Need to change if add more
J_params = [1 2];
p_params = [4 5]; 

for icond = 1:length(condVec)
    cond = condVec{icond};
    %select Jbar, p_high, and data based on TMS condition
    this_theta = Theta([J_params(icond),3,p_params(icond)]); 
    this_data = data.(cond);
    
    nLL = nLL + calc_nLL('flexible',this_theta,this_data,exppriorityVec);
end
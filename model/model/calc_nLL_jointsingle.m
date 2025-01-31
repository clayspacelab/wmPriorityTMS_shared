function nLL = calc_nLL_jointsingle(Theta,data,exppriorityVec,fixparams,condVec)
%fits the basic VP model separately but jointly to all TMS/priority conditions

nLL = 0;

for icond = 1:length(condVec)
    cond = condVec{icond};
    
    %select data based on TMS condition and priority
    for ip = 1:length(exppriorityVec)
        this_data = data.(cond){:,ip};
        this_thetastart = (icond-1)*4+(ip-1)*2+1; %sequence thetas in order (ugly)
        this_theta = Theta(this_thetastart:(this_thetastart+1));
        
        nLL = nLL + calc_nLL_single('single',this_theta,this_data);
    end
end
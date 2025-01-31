function nLL = calc_nLL_base(Theta,data,exppriorityVec,fixparams,condVec)
%assumes no effect of TMS (e.g. just fit the basic 3 parameter priority
%model to all data)

nLL = 0;

for icond = 1:length(condVec)
    cond = condVec{icond};
    %select data based on TMS condition
    this_data = data.(cond);
    
    nLL = nLL + calc_nLL('flexible',Theta,this_data,exppriorityVec);
end
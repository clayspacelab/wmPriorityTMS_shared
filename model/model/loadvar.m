function varargout = loadvar(varargin)
nvars = ceil(length(varargin)/2);
varargout = cell(1,nvars);
for ivar = 1:nvars
    var = varargin{2*ivar-1};
    switch var
        case 'JVec'
            Jbar = varargin{2*ivar}{1};
            tau = varargin{2*ivar}{2};
            %{
            nJSamp = 500;
            
            JVec = linspace(1e-10,10*Jbar,nJSamp);
            %}
            
            nJSamp = 1000;
            %prange = linspace(1e-4,1-1e-4,nJSamp);
            
           
            %sample probabilities in logit space, 
            %then use logistic and inverse gamma to get Js(best method at
            %adequately sampling Js that I've found, as it produces fewest
            %discretizations that don't integrate to 1)
            prangel = linspace(log(1e-4./(1-1e-4)),log((1-1e-4)./1e-4),nJSamp);
            prange = 1./(1+exp(-prangel));
            JVec = gaminv(prange,Jbar/tau,tau);
            
            %{
            %sampling more in the middle probabilities inuitively seems better but doesn't work for
            %having pdf that integrates to 1
            prange = linspace(1e-4,1-1e-4,nJSamp);
            JVec = gaminv(prange,Jbar/tau,tau);
            %}

            %{
            %sample Js using inverse gamma for endpoints, skipping
            %probability step. Also didn't work well
            %JVec = linspace(gaminv(1e-4,Jbar/tau,tau),gaminv(1-1e-4,Jbar/tau,tau),nJSamp);
            %prange = gamcdf(prange2,Jbar/tau,tau);
            %}
            
            
            %xgampdf = gampdf(JVec,Jbar/tau,tau);
            %trapz(JVec,xgampdf)
            
%             nJSamp = 100;
%             xmin = Jbar;
%             xmax = Jbar;
% 
%                 % lower bound
%                 while gampdf(xmin,Jbar/tau,tau) > 1e-4
%                     if xmin <= 1
%                         xmin = 1e-10;
%                         break
%                     else
%                         xmin = xmin - 1;
%                     end
%                 end
%                 
%                 increment = 1;
%                 % upper bound
%                 while gampdf(xmax,Jbar/tau,tau) > 1e-4
%                     xmax = xmax + increment;
%                 end
%                 % decrease by smaller increments to make it more fine
%                 % grained
%                 increment = 0.1*increment;
%                 while gampdf(xmax,Jbar/tau,tau) < 1e-4
%                     xmax = xmax - increment;
%                     if xmax <= 0.1
%                         xmax = 0.05;
%                         break
%                     end
%                 end
% %             end
%             JVec = linspace(xmin,xmax,nJSamp);
% %             JVec = linspace(1e-5,20,nJSamp); % ASPEN: make sure range is reasonable for parameter range
            varargout{ivar} = JVec;
        case 'rVec'
%             nRs = varargin{2*ivar};
            % radius stuff
            nRs = 500;
            rVec = linspace(0,10,nRs);
            varargout{ivar} = rVec;
    end
end
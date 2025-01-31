function nLL = calc_nLL(model,Theta,data,exppriorityVec,fixparams)
%CALC_NLL calculates negative log-likelihood of parameters given data and
%model.
% 
%   NLL = CALC_NLL(MODEL, THETA, DATA, EXPPRIORITYVEC) calculates the 
%     negative likelihood of DATA given MODEL and THETA, for experimental
%     probe probabilities EXPPRIORITYVEC.
% 
%   NLL = CALC_NLL(MODEL, THETA, DATA, EXPPRIORITYVEC, FIXPARAMS) 
%     calculates the negative likelihood, where THETA are free parameters 
%     and FIXPARAMS indicates which parameters are fixed and what value 
%     they are fixed to. 
%
%   ================= INPUT VARIABLES ======================
% 
%   MODEL / THETA: 
%         'max_points', MAXIMIZING POINTS / [Jbar_total tau (alpha beta)]
%         'flexible', FLEXIBLE / [Jbar_total tau (alpha beta) p_high p_med]
%         'proportional', PROPORTIONAL / [Jbar_total tau (alpha beta)]
%         'min_error', MINIMIZING ERROR / [Jbar_total tau (alpha beta) gamma]
% 
%       parameter descriptions: 
%           JBAR_TOTAL: mean total amount of resources across priorities
%           TAU: second parameter of gamma noise distribution
%           ALPHA: risk preferences for post-decision wager
%           BETA: inverse noise temperature on post-decision wager
%           GAMMA: exponent for loss in Minimizing Error model
%           P_HIGH: proportion allocated to high-priority stimulus
%           P_MED: proportion allocated to medium-priority stimulus
%
%   DATA: cell of length nPriorities. the ith element of DATA should be
%     all data corresponding to EXPPRIORITYVEC(i) condition. the first
%     column should contain the magnitude of errors and the second column,
%     if available, should contain the corresponding circle wager radius
%     size. 
% 
%   FIXPARAMS: (optional). 2 x (number of fixed parameters) matrix. fixed 
%     parameters, such that the first row corresponds to the index and 
%     second row corresponds to the value of the fixed parameter. 
%
%   ================= OUTPUT VARIABLES ================
% 
%   NLL: negative log-likelihood

% ---------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% ---------------------

if nargin < 5; fixparams = []; end

expnumber = size(data{1},2); % experiment number

% exponentiating appropriate parameters
logflag = loadconstraints(model,exppriorityVec,expnumber-1);

% if there are fixed parameters
if ~isempty(fixparams)
    logflag(fixparams(1,:)) = 0;
    
    nParams = length(Theta) + size(fixparams,2);
    nonfixedparamidx = 1:nParams;
    nonfixedparamidx(fixparams(1,:)) = [];
    
    temptheta = nan(1,nParams);
    temptheta(nonfixedparamidx) = Theta;
    temptheta(fixparams(1,:)) = fixparams(2,:);
    
    Theta = temptheta;
end

Theta(logflag) = exp(Theta(logflag));

Jbar_total = Theta(1);
tau = Theta(2);
if (expnumber == 2); alpha = Theta(3); beta = Theta(4); end

%set optional noise parameter, which can be used to model convolving the
%memory noise with a sensorimotor noise component
if strcmp(model,'flexible_anoise') && expnumber == 1 %just to play nice
    anoise = Theta(3);
else
    anoise = 0;
end

% Theta

% obtain vector of resource allocated
nPriorities = length(exppriorityVec);
switch model
    case 'max_points'   % maximizing points (exp 2 only)
        pVec = calc_pVec_maxpoints(Theta,exppriorityVec);
    case {'flexible','flexible_anoise'}     % flexible
        pp = Theta(end-(nPriorities-2):end);
        pVec = [pp 1-sum(pp)];
    case 'proportional' % proportional
        pVec = exppriorityVec;
    case 'min_error'    % minimizing error^\gamma
        pVec = calc_pVec_minerror(Theta,exppriorityVec);
        if (roundn(sum(pVec),-3) ~=1)
            fprintf('For theta = [%.3f %.3f %.3f], pVec does not add to one: [%.3f %.3f] \n',...
                Theta(1), Theta(2), Theta(3), pVec(1),pVec(2));
        end
    case 'flex_jbars'
        Jbar_total = sum(Theta(1:2));
        tau = Theta(3);
        pVec = Theta(1:2)./Jbar_total;
end


% loading vector of circle wager radii
[rVec] = loadvar('rVec'); % size: (1 x 500)
rVec = rVec(:); % size: (500 x 1)

if any(isinf(pVec))
    nLL = Inf;
else
    nLL = 0;
    for ipriority = 1:nPriorities
        
        idx_del = find(sum(isnan(data{ipriority}),2));
        data{ipriority}(idx_del,:) = [];
        Jbar = Jbar_total*pVec(ipriority); % Jbar for current priority condition
        
        % clear varaibles used in previous priorities (necessary for code to run)
        clear idx1 idx2 data_r_reshaped
        
        % get subject data
        data_distance = data{ipriority}(:,1);
        nTrials = length(data_distance);
        
        % p(J|Jbar,tau)
        [JVec] = loadvar('JVec',{Jbar,tau}); % values of J
        nJs = length(JVec);
        Jpdf = gampdf(JVec,Jbar/tau,tau); % probability density over Js
        Jpdf_int = trapz(JVec,Jpdf);
        if abs(Jpdf_int-1) > 1e-2
            warning('pdf does not integrate to 1: %f (J=%f,tau=%f)!',Jpdf_int,Jbar,tau)
        end
        %assert(trapz(JVec,Jpdf)-1 < 1e-2,'pdf does not integrate to 1: %f (J=%f,tau=%f)!',trapz(JVec,Jpdf),Jbar,tau)
        Jpdf = Jpdf./Jpdf_int; % normalize slop
        %figure();plot(JVec,Jpdf)
        %{
        Jpdf2 = Jpdf;
        Jpdf = Jpdf./sum(Jpdf); % normalize
        trapz(JVec,Jpdf2)
        assert(iscloseall(trapz(JVec,Jpdf2),1),'pdf does not integrate to 1: %f!',trapz(JVec,Jpdf2))
        %}
        
        % p(Shat|S,J)
        %{
        Sigma = zeros(1,2,nJs*nTrials);
        Sigma(1,:,:) = sort(repmat(sqrt(1./JVec(:)),nTrials,2),'descend')'; % SDs of diagonal matrix. sigmas in descending order --> J in ascending order
        p_Shat = mvnpdf(repmat([data_distance(:) zeros(nTrials,1)],nJs,1),0,Sigma);
        p_Shat = reshape(p_Shat,nTrials,nJs)'; % nJs x nTrials
        %}
        %(p_err|J)
        p_Shat = raylpdf(repmat(data_distance',nJs,1),...
            repmat(sqrt(1./JVec) + anoise,nTrials,1)');
        p_Shat(p_Shat == 0) = 1e-10; % set to arbitrarily small value if zero
        
        % ====== Exp 2: with disc size data ======
        if (expnumber == 2)
            data_r = data{ipriority}(:,2);
            
            % p(rVec|J,beta) (a range of r to get entire probability dist)
            pdf_r = calc_pdf_r(beta, JVec, alpha); % size: (nrs x nJs)
            
            xdiff = diff(rVec(1:2));
            p_r = nan(nJs,nTrials);
            for iJ = 1:nJs
                pdff = pdf_r(:,iJ);
                for itrial = 1:nTrials
                    r = data_r(itrial);
                    idx1 = find((rVec- r) <= 0,1,'last');
                    idx2 = find((rVec- r) > 0,1,'first');
                    
                    p_r(iJ,itrial) = (pdff(idx2)-pdff(idx1))/xdiff.*(r-rVec(idx1)) + pdff(idx1);
                end
            end
        end
        
        if (expnumber == 2) % if there is wager data
            % \int p(Shat|S,J) p(r|J) p(J) dJ
            pTrials = sum(bsxfun(@times,p_Shat.*p_r,Jpdf')); % 1 x nTrials
        else
            % \int p(Shat|S,J) p(J) dJ
            %pTrials = sum(bsxfun(@times,p_Shat,Jpdf')); % 1 x nTrials
            pTrials = trapz(JVec',p_Shat.*Jpdf'); % 1 x nTrials %marginalize over Js
        end
        %pTrials3 = trapz(JVec',p_Shat.*Jpdf2');
        %assert(iscloseall(pTrials,pTrials3),'trial prs do not match')
        %figure();scatter(pTrials3,pTrials);hold on;fplot(@(x) x,[min(pTrials3) max(pTrials3)])
        nLL = nLL - sum(log(pTrials));
        
    end
end
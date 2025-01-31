function nLL = calc_nLL_single(model,Theta,data)
% calc_nLL_single calculates VP model for a single dataset (i.e, no
% priority conditions, etc.)
%   ================= OUTPUT VARIABLES ================
% 
%   NLL: negative log-likelihood

% exponentiating appropriate parameters
logflag = loadconstraints(model,[],0);
Theta(logflag) = exp(Theta(logflag));

Jbar = Theta(1);
tau = Theta(2);

%idx_del = find(isnan(data));
%data(idx_del,:) = [];
data(isnan(data),:) = [];

% get length of data
nTrials = length(data);

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
p_Shat = mvnpdf(repmat([data(:) zeros(nTrials,1)],nJs,1),0,Sigma);
p_Shat = reshape(p_Shat,nTrials,nJs)'; % nJs x nTrials
%}
%(p_err|J)
p_Shat = raylpdf(repmat(data',nJs,1),repmat(sqrt(1./JVec),nTrials,1)');
p_Shat(p_Shat == 0) = 1e-10; % set to arbitrarily small value if zero

% \int p(Shat|S,J) p(J) dJ
%pTrials = sum(bsxfun(@times,p_Shat,Jpdf')); % 1 x nTrials
pTrials = trapz(JVec',p_Shat.*Jpdf'); % 1 x nTrials %marginalize over Js

%pTrials3 = trapz(JVec',p_Shat.*Jpdf2');
%assert(iscloseall(pTrials,pTrials3),'trial prs do not match')
%figure();scatter(pTrials3,pTrials);hold on;fplot(@(x) x,[min(pTrials3) max(pTrials3)])
nLL = -sum(log(pTrials));
        
end
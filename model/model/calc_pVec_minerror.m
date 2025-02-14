function [pVec, fval] = calc_pVec_minerror(Theta,exppriorityVec)
%CALC_PVEC_MINERROR calculates the proportion allocated to each priority 
%condition that minimizes loss
%
%   =========== INPUT VARIABLES =========
% 
%   THETA: [Jbar_total, tau, [alpha, beta,] gamma]. 
% 
%   EXPPRIORITYVEC: 1 x 3 vector of experimental priority.
%   sum(exppriorityVec) = 1

nPriorities = length(exppriorityVec)-1;

% function for expected error
calc_E_error = @(x) calc_expectederror_analytical(Theta,[x 1-sum(x)],exppriorityVec);

% parameters for optimization
Aeq = [];%ones(1,nPriorities); 
beq = [];%1; 
A = ones(1,nPriorities); %diag(ones(1,nPriorities).*(-Theta(1)/Theta(2))); 
b = 1; %-Theta(end)/2.*ones(nPriorities,1); 
nonlcon = deal([]); 
options = optimset('Display','none'); 
lb = 1e-3.*ones(1,nPriorities);
ub = ones(1,nPriorities);
nStartVals = 10;  % tried with different parameters and lowest value showed up 3,5,7,8,10,10 of 10. 

% lower and upper bounds for starting pts for pVec
lbb = ones(1,nPriorities).*(Theta(end)*Theta(2)/Theta(1)/2); % ASPEN WHY IS THIS LOWER BOUND
ubb = 1-sum(lbb);

if lbb(1) > ubb
    u = ubb;
    ubb = lbb(1);
    lbb = ones(1,nPriorities-1).*u;
end


if (nPriorities == 1)
    x0 = rand(nStartVals,1);
else
    x0 = [];
    constantt = 0;
    while size(x0,1) < nStartVals
        x0 = lhs(nStartVals+constantt,nPriorities-1,lbb,ubb,[],1e3);
        x0 = [x0 1-sum(x0,2)];      % all priorities sum to 1
        idx = (x0(:,end) < lbb(1)); % delete x0's when they aren't within the bounds
        x0(logical(idx),:) =[];
        constantt = constantt + sum(idx);
    end
end
% optimizing
pVec = nan(nStartVals,nPriorities);
nEU = nan(1,nStartVals);

tic;
for istartval = 1:nStartVals
    [pVec(istartval,:), nEU(istartval)] = fmincon(calc_E_error,x0(istartval,:),A,b,Aeq,beq,lb,ub,nonlcon,options);
    if toc > 30 % seconds
        fval = [];
        fprintf('estimating minimizing error resource allocation took longer than 30 seconds \n')
        return
    end
end

pVec(nEU < 0,:) = [];
nEU(nEU < 0) = [];

fval = min(nEU);
pVec = pVec(nEU == fval,:);

if isempty(fval)
    fval = Inf;
    pVec = Inf*ones(1,nPriorities);
else
    pVec = pVec(1,:);  % in case multiple entries have the nEU == min(nEU)
end

pVec = [pVec 1-sum(pVec)];

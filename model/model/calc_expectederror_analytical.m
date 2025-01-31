function expectederror = calc_expectederror_analytical(Theta,allocatedpriorityVec,exppriorityVec)
%CALC_EXPECTEDERROR_ANALYTICAL analytically computes the expected cost
%
%   CALC_EXPECTEDERROR_ANALYTICAL(THETA,ALLOCATEDPRIORITYVEC) calculates the
%     expected cost (euclidean error ^ gamma) for a given parameter vector THETA
%     and resource allocation across consitions ALLOCATEDPRIORITYVEC.
% 
%   ===== INPUT VARIABLES =====
%   THETA: Jbar_total, tau, [alpha, beta,] gamma. (alpha and beta are only in
%     experiment 2)
%   ALLOCATEDPRIORITYVEC: row vector of allocated priority.
%     sum(allocatedpriorityVec) = 1
%   EXPPRIORITYVEC: row vector of experimental priority.
%     sum(exppriorityVec) = 1

% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% -----------------------


% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
gamma = Theta(end);

nPriorities = length(allocatedpriorityVec);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
%            = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

expectederror = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);

    k = Jbar/tau;
    
    if any([gamma/2+1 k-(gamma/2) k] <= 0) % this constraint must be met for expected error to be analytical
        % a note to self and to potential readers: i have tried to
        % implement a numerical solution here, but the estimates are a bit
        % off. i am not sure if it is better to have an estimate and be off
        % than no estimate at all, but right now it is the latter. see
        % CALC_EXPECTEDERROR_NUMERICAL.M and
        % CALC_EXPECTEDERROR_HALFNUMERICAL.M for the current attempts at
        % solutions for this issue. 
        
        % september 16, 2022: In grace's study, running into a problem
        % where very high gammas are put into this function and break,
        % trying out sticking numerical solution here, recognizing it may
        % be incorrect. 
        bleh = calc_expectederror_numerical(Theta,allocatedpriorityVec(ipriority),exppriorityVec(ipriority));
        expectederror = expectederror + bleh;
%         Theta
%         allocatedpriorityVec
%         sprintf('estimated numerical cost for allocation %.2f is %02.4f \n',allocatedpriorityVec(ipriority),bleh/exppriorityVec(ipriority))
%         expectederror = Inf
%         return
    else
        bleh = exp(gammaln(gamma/2 + 1) + gammaln(k-(gamma/2)) - gammaln(k)) .* (2/tau)^(gamma/2);
        expectederror = expectederror + exppriorityVec(ipriority).*bleh;
%         sprintf('estimated analytical cost for allocation %.2f is %02.4f \n',allocatedpriorityVec(ipriority),bleh)
    end    
end


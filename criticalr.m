function [crit_r] = criticalr(N,alpha)

if nargin==1
    alpha=0.05;
end;

DF= N-2; % degrees of freedom
crit_t = nctinv(alpha/2,DF,0); % quantiles for the noncentral student t distribution
crit_r = sqrt((crit_t^2)/((crit_t^2)+DF)); % critical r value
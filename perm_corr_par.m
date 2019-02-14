function [real_res, p_real, s_real, h_real, dummy_res, dummy_hvs] = perm_corr_par(data, covd, type, alpha, nPerm)
%% [real_res, dummy_res] = corr_perm(data, covd, type, nPerm)
% Giulio Bernardi [giulioberna@gmail.com], 2017.11.20

disp('Calculation of real and null data values...');
[real_res, p_real]=corr(data,covd,'type',type,'rows','pairwise');
h_real=p_real<alpha;
s_real=h_real; s_real(real_res<0 & h_real==1)=-1;
dummy_res=NaN(size(real_res,1),nPerm);
dummy_hvs=NaN(size(real_res,1),nPerm);
for i=1:nPerm  
    if mod(i, 100)==0; disp(['  ',num2str(i),' out of ',num2str(nPerm),' permutations completed']); end
    rand_covd=covd(randperm(length(covd)));
    [dummy_res(:,i), pvals]=corr(data,rand_covd,'type',type,'rows','pairwise');
    dummy_hvs(:,i)=(pvals<alpha).*(sign(dummy_res(:,i))); clear pvals;
end; clear i;

end
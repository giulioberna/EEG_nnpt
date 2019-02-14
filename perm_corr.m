function [real_res, dummy_res] = perm_corr(data, covd, type, nPerm)
%% [real_res, dummy_res] = corr_perm(data, covd, type, nPerm)
% Giulio Bernardi [giulioberna@gmail.com], 2017.11.19

disp('Calculation of real and null data values...');
real_res=corr(data,covd,'type',type,'rows','pairwise');
dummy_res=NaN(size(real_res,1),nPerm);
for i=1:nPerm  
    if mod(i, 100)==0; disp(['  ',num2str(i),' out of ',num2str(nPerm),' permutations completed']); end
    rand_covd=covd(randperm(length(covd)));
    dummy_res(:,i)=corr(data,rand_covd,'type',type,'rows','pairwise');
end; clear i;

end
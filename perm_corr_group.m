function [real_res, dummy_res] = perm_corr_group(data_sub, covd_sub, type, nPerm)
%% [real_res, dummy_res] = perm_corr_group_par(data, covd, type, nPerm)
% data_sub and covd_sub are cells; data_sub is organized as observation x channel
% Giulio Bernardi [giulioberna@gmail.com], 2017.12.21

disp('Calculation of real and null data values...');

nsubs=length(data_sub);

rho_real=[]; 
rho_dummy=[];

for ss=1:nsubs
    
    disp(['    >> Subject ',num2str(ss)]);

    data=data_sub{ss};
    covd=covd_sub{ss};
    
    % Real Results
    if strcmp(type,'Spearman')
        [rho_real(:,ss)]=corr(data,covd,'rows','complete','type','Spearman');
    elseif strcmp(type,'Pearson')
        [temp_rho]=corr(data,covd,'rows','complete','type','Pearson');
        for ch=1:length(temp_rho)
            if abs(temp_rho(ch))~=1
                rho_real(ch,ss)=fisherz(temp_rho(ch)); 
            else
                rho_real(ch,ss)=NaN;
            end;
        end; clear temp_rho ch;
    end;
    
    % Dummy Results
    for i=1:nPerm  
        if mod(i, 100)==0; disp(['  ',num2str(i),' out of ',num2str(nPerm),' permutations completed']); end
        rand_covd=covd(randperm(length(covd)));
        if strcmp(type,'Spearman')
            [rho_dummy(:,ss,i)]=corr(data,rand_covd,'rows','complete','type','Spearman');
        elseif strcmp(type,'Pearson')
            [temp_rho]=corr(data,rand_covd,'rows','complete','type','Pearson');
            for ch=1:length(temp_rho)
                if abs(temp_rho(ch))~=1
                    rho_dummy(ch,ss,i)=fisherz(temp_rho(ch,i)); 
                else
                    rho_dummy(ch,ss,i)=NaN;
                end;
            end; clear temp_rho ch;
        end;  
    end; clear i;
    
end; clear ss;

real_res=nanmean(rho_real,2); 
dummy_res=squeeze(nanmean(rho_dummy,2));

end
  

%%
function[z]=fisherz(r)
% Fisher's Z-transform.

r=r(:);
z=.5.*log((1+r)./(1-r));

end
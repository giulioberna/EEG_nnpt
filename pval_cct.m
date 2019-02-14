function [t_clust, h_real_cc] = pval_cct(h_real,s_real,dummy_res,NBRS,tail,alpha,calpha)
%% [t_clust, h_real_cc] = pval_cct(h_real,s_real,dummy_res,NBRS,tail,alpha,calpha)
% Giulio Bernardi [giulioberna@gmail.com], 2017.11.19
% 
% This function requires the surfing_clusterize function from the
% surfing-master toolbox.

disp('Calculation of cluster-based correction...');
h_real_cc=zeros(size(dummy_res,1),1); 
dummy_pvs=NaN(size(dummy_res)); 
dummy_sde=NaN(size(dummy_res));
dummy_mcl=NaN(2,size(dummy_res,2)); 
for i=1:size(dummy_pvs,2)
    if mod(i, 100)==0; disp(['  ',num2str(i),' out of ',num2str(size(dummy_res,2)),' iterations completed']); end
    for d=1:size(dummy_pvs,1)
        temp_right = (length(find(dummy_res(d,:)>=dummy_res(d,i)))-1)/(size(dummy_res,2)-1);
        temp_left  = (length(find(dummy_res(d,:)<=dummy_res(d,i)))-1)/(size(dummy_res,2)-1);
        if strcmpi(tail,'right') || (strcmpi(tail,'both') && dummy_res(d,i) >= median(dummy_res(d,:)))
            dummy_pvs(d,i) = temp_right;
            dummy_sde(d,i) = 1;
        elseif strcmpi(tail,'left') || (strcmpi(tail,'both') && dummy_res(d,i) <= median(dummy_res(d,:)))
            dummy_pvs(d,i) = temp_left;
            dummy_sde(d,i) = -1;
        end; 
    end; clear d; 
    temp_dummy=(dummy_pvs(:,i)<alpha).*(dummy_sde(:,i)==1);
    temp_clust=surfing_clusterize(temp_dummy,NBRS);
    temp_csize=sort(cell2mat(cellfun(@(x)(length(x)),temp_clust,'uniformoutput',false)));
    if isempty(temp_csize); dummy_mcl(1,i)=0; else, dummy_mcl(1,i)=max(temp_csize); end;
    temp_dummy=(dummy_pvs(:,i)<alpha).*(dummy_sde(:,i)==-1);
    temp_clust=surfing_clusterize(temp_dummy,NBRS);
    temp_csize=sort(cell2mat(cellfun(@(x)(length(x)),temp_clust,'uniformoutput',false)));
    if isempty(temp_csize); dummy_mcl(2,i)=0; else, dummy_mcl(2,i)=max(temp_csize); end;
    clear temp_dummy temp_clust temp_csize temp_right temp_left;
end; clear i;

idx=[];
switch tail 
    case 'right'
        t_clust(1)=prctile(dummy_mcl(1,:),(1-calpha)*100);
        temp_data=h_real.*(s_real==1);
        temp_clust=surfing_clusterize(temp_data,NBRS);
        temp_csize=cell2mat(cellfun(@(x)(length(x)),temp_clust,'uniformoutput',false));
        idx=[idx;cat(1,temp_clust{temp_csize>t_clust(1)},1)]; clear temp_data temp_clust temp_csize;
    case 'left' 
        t_clust(2)=prctile(dummy_mcl(2,:),(1-calpha)*100);
        temp_data=h_real.*(s_real==-1);
        temp_clust=surfing_clusterize(temp_data,NBRS);
        temp_csize=cell2mat(cellfun(@(x)(length(x)),temp_clust,'uniformoutput',false));
        idx=[idx;cat(1,temp_clust{temp_csize>t_clust(2)})]; clear temp_data temp_clust temp_csize;
    case 'both' 
        t_clust(1)=prctile(dummy_mcl(1,:),(1-(calpha/2))*100);
        temp_data=h_real.*(s_real==1);
        temp_clust=surfing_clusterize(temp_data,NBRS);
        temp_csize=cell2mat(cellfun(@(x)(length(x)),temp_clust,'uniformoutput',false));
        idx=[idx;cat(1,temp_clust{temp_csize>t_clust(1)})]; clear temp_data temp_clust temp_csize;
        t_clust(2)=prctile(dummy_mcl(2,:),(1-(calpha/2))*100);
        temp_data=h_real.*(s_real==-1);
        temp_clust=surfing_clusterize(temp_data,NBRS);
        temp_csize=cell2mat(cellfun(@(x)(length(x)),temp_clust,'uniformoutput',false));
        idx=[idx;cat(1,temp_clust{temp_csize>t_clust(2)})]; clear temp_data temp_clust temp_csize;
end;
h_real_cc(idx)=1; % cluster-based correction

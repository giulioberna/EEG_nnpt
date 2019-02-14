function [nppt_corr_results, dummy_res] = nppt_corr_group(data, covd, type, NBRS, nPerm, tail, alpha, calpha)
%% [nppt_corr_results] = nppt_cor_parcc(data, covd, type, NBRS, nPerm, tail, alpha, calpha)
% Permutation-based correlation analysis incorporating correction for multiple comparisons
% For group-level correction on correlations performed at single-subject level
% Giulio Bernardi [giulioberna@gmail.com], 2017.12.21
    
    %% Data quality checks
    
    if size(data,1) ~= size(covd,1)
        disp('ERROR: Data and covariate should have the same number of observations!')
    return;
    end;
    
    disp(['Permutation Test with base p=',num2str(alpha),' and corrected p=', num2str(calpha)]);
    disp(['   Correlation Analysis with tail: ',tail]);
    disp(['   Number of Permutations: ',num2str(nPerm)]);
    
    tic;
    
    nch=size(data{1},2);
    
    %% Initialization
    nppt_corr_results.options.type=type;               % used options (type)
    nppt_corr_results.options.nPerm=nPerm;             % used options (nPerm)
    nppt_corr_results.options.tail=tail;               % used options (tail)
    nppt_corr_results.options.alpha=alpha;             % used options (alpha)
    nppt_corr_results.options.calpha=calpha;           % used options (calpha)
    nppt_corr_results.v_real=NaN(nch,1);               % results of the correlation
    nppt_corr_results.p_real=NaN(nch,1);               % uncorrected p-values
    nppt_corr_results.s_real=NaN(nch,1);               % side of the comparison
    nppt_corr_results.z_real=NaN(nch,1);               % z-score
    nppt_corr_results.h_real=NaN(nch,1);               % hypothesis test uncorrected
    nppt_corr_results.h_real_cc=zeros(nch,1);          % cluster-based correction
    nppt_corr_results.h_real_pc=zeros(nch,1);          % p-based correction
    nppt_corr_results.t_pcorr=NaN(1,2);                % threshold for p-correction
    nppt_corr_results.t_clust=NaN(1,2);                % threshold for c-correction

    %% Calculation of real and null data values
    [real_res, dummy_res] = perm_corr_group(data, covd, type, nPerm);
    nppt_corr_results.v_real=real_res; 
    
    %% Calculation of uncorrected p-values for real data
    [nppt_corr_results.p_real, nppt_corr_results.s_real, nppt_corr_results.z_real] = pval_unt(real_res,dummy_res,tail);
    nppt_corr_results.h_real=nppt_corr_results.p_real<alpha; 
    
     %% P-based correction
    [nppt_corr_results.t_pcorr, nppt_corr_results.h_real_pc] = pval_pct(real_res,dummy_res,tail,calpha);
          
    %% Cluster-based correction
    [nppt_corr_results.t_clust, nppt_corr_results.h_real_cc] = pval_cct(nppt_corr_results.h_real,nppt_corr_results.s_real,dummy_res,NBRS,tail,alpha,calpha);
    
    %% Send output on command window
    disp('Operation completed...');
    disp(['  ',num2str(sum(nppt_corr_results.h_real)),' values significant at uncorrected level [',num2str(round(sum(nppt_corr_results.h_real)./length(nppt_corr_results.h_real)*100)),'%]']);
    disp(['  ',num2str(sum(nppt_corr_results.h_real_pc)),' values significant at p-based corrected level [',num2str(round(sum(nppt_corr_results.h_real_pc)./length(nppt_corr_results.h_real_pc)*100)),'%]']);
    disp(['  ',num2str(sum(nppt_corr_results.h_real_cc)),' values significant at cluster-based corrected level [',num2str(round(sum(nppt_corr_results.h_real_cc)./length(nppt_corr_results.h_real_cc)*100)),'%]']);
    
    toc; % Display Elapsed Time
    
end % End of Function
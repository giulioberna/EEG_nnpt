function [nppt_test_results, dummy_res] = nppt_test(data1, data2, type, NBRS, nPerm, tail, alpha, calpha)
%% [nppt_test_results] = nppt_test(data1, data2, type, NBRS, nPerm, tail, alpha, calpha)
% Permutation-based contrast analysis incorporating correction for multiple comparisons
% Giulio Bernardi [giulioberna@gmail.com], 2018.09.07
    
    %% Data quality checks
    
    if size(data1,1) ~= size(data2,1)
        disp('ERROR: The two groups should have the same number of variables!')
    return;
    end;
    
    disp(['Permutation Test with base p=',num2str(alpha),' and corrected p=', num2str(calpha)]);
    disp(['   ',type,' Group Comparison with tail: ',tail]);
    disp(['   Number of Permutations: ',num2str(nPerm)]);
    
    tic;
    
    %% Initialization
    nppt_test_results.options.type=type;                % used options (type)
    nppt_test_results.options.nPerm=nPerm;              % used options (nPerm)
    nppt_test_results.options.tail=tail;                % used options (tail)
    nppt_test_results.options.alpha=alpha;              % used options (alpha)
    nppt_test_results.options.calpha=calpha;            % used options (calpha)
    nppt_test_results.v_real=NaN(size(data1,2),1);      % results of the correlation
    nppt_test_results.p_real=NaN(size(data1,2),1);      % uncorrected p-values
    nppt_test_results.s_real=NaN(size(data1,2),1);      % side of the comparison
    nppt_corr_results.z_real=NaN(size(data1,2),1);      % z-score
    nppt_test_results.h_real=NaN(size(data1,2),1);      % hypothesis test uncorrected
    nppt_test_results.h_real_cc=zeros(size(data1,2),1); % cluster-based correction
    nppt_test_results.h_real_pc=zeros(size(data1,2),1); % p-based correction
    nppt_test_results.t_pcorr=NaN(1,2);                 % threshold for p-correction
    nppt_test_results.t_clust=NaN(1,2);                 % threshold for c-correction

    %% Calculation of real and null data values
    [real_res, dummy_res] = perm_test(data1, data2, type, nPerm);
    nppt_test_results.v_real=real_res;
        
    %% Calculation of uncorrected p-values for real data
    [nppt_test_results.p_real, nppt_test_results.s_real, nppt_test_results.z_real] = pval_unt(real_res,dummy_res,tail);
    nppt_test_results.h_real=nppt_test_results.p_real<alpha; 
    
    %% P-based correction
    [nppt_test_results.t_pcorr, nppt_test_results.h_real_pc] = pval_pct(real_res,dummy_res,tail,calpha);
    
    %% Cluster-based correction
    [nppt_test_results.t_clust, nppt_test_results.h_real_cc] = pval_cct(nppt_test_results.h_real,nppt_test_results.s_real,dummy_res,NBRS,tail,alpha,calpha);
    
    %% Send output on command window
    disp('Operation completed...');
    disp(['  ',num2str(sum(nppt_test_results.h_real)),' values significant at uncorrected level [',num2str(round(sum(nppt_test_results.h_real)./length(nppt_test_results.h_real)*100)),'%]']);
    disp(['  ',num2str(sum(nppt_test_results.h_real_pc)),' values significant at p-based corrected level [',num2str(round(sum(nppt_test_results.h_real_pc)./length(nppt_test_results.h_real_pc)*100)),'%]']);
    disp(['  ',num2str(sum(nppt_test_results.h_real_cc)),' values significant at cluster-based corrected level [',num2str(round(sum(nppt_test_results.h_real_cc)./length(nppt_test_results.h_real_cc)*100)),'%]']);
    
    toc; % Display Elapsed Time
    
end % End of Function
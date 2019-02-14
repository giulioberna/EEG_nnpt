function [real_res, dummy_res] = perm_slope(data, dcov, nPerm)
%% [real_res, dummy_res] = test_perm(data1, data2, type, nPerm)
% Giulio Bernardi [giulioberna@gmail.com], 2017.11.21

disp('Calculation of real and null data values...');

sub_res=NaN(size(data,1),size(data,2));
for s=1:size(data,2)
        temp_data=squeeze(data(:,s,:)); 
        temp_data=temp_data(:,~isnan(temp_data(1,:)));
        temp_dcov=dcov(s,~isnan(dcov(s,:)));
        if length(temp_dcov)>1
            for d=1:size(temp_data,1) 
                [p] = calcslope(temp_dcov,temp_data(d,:),2,0);
                sub_res(d,s)=p(2);
            end; clear d;
        end; clear d;
end; clear s;
real_res=nanmedian(sub_res,2);

dummy_res=NaN(size(real_res,1),nPerm);
for i=1:nPerm
    if mod(i, 100)==0; disp(['  ',num2str(i),' out of ',num2str(nPerm),' permutations completed']); end
    sub_res=NaN(size(data,1),size(data,2));
    for s=1:size(data,2)
            temp_data=squeeze(data(:,s,:)); 
            temp_data=temp_data(:,~isnan(temp_data(1,:)));
            temp_dcov=dcov(s,~isnan(dcov(s,:)));
            temp_dcov=temp_dcov(randperm(length(temp_dcov)));
            if length(temp_dcov)>1
                for d=1:size(temp_data,1) 
                    [p] = calcslope(temp_dcov,temp_data(d,:),2,0);
                    sub_res(d,s)=p(2);
                end; clear d;
            end; clear d;
    end; clear s;
    dummy_res(:,i)=nanmedian(sub_res,2);
end; clear i;

end
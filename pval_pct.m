function [t_pcorr, h_real_pc] = pval_pct(real_res,dummy_res,tail,calpha)
%% [t_pcorr, h_real_pc] = pval_pct(real_res,dummy_res,tail,calpha)
% Giulio Bernardi [giulioberna@gmail.com], 2017.11.19

disp('Calculation of p-based correction...');
t_pcorr=NaN(1,2); h_real_pc=zeros(size(real_res));
dummy_mpv=[max(dummy_res,[],1);min(dummy_res,[],1)];
switch tail 
    case 'right'
        t_pcorr(1)=prctile(dummy_mpv(1,:),(1-calpha)*100);
        h_real_pc(real_res>t_pcorr(1))=1; % p-based correction
    case 'left'
        t_pcorr(2)=prctile(dummy_mpv(2,:),calpha*100);
        h_real_pc(real_res<t_pcorr(2))=1; % p-based correction
    case 'both' 
        t_pcorr(1)=prctile(dummy_mpv(1,:),(1-(calpha/2))*100);  
        h_real_pc(real_res>t_pcorr(1))=1; % p-based correction
        t_pcorr(2)=prctile(dummy_mpv(2,:),(calpha/2)*100);
        h_real_pc(real_res<t_pcorr(2))=1; % p-based correction
end;

end
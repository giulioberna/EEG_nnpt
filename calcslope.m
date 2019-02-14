function [p] = calcslope(x,y,type,fig)

if type==1 % Robust Fit
    p = robustfit(x,y);
elseif type==2 % Polyfit
    p = polyfit(x,y,1);
    p = fliplr(p);
end;

if fig==1
    figure; scatter(x,y,'fill'); hold on; plot(x,p(2)*x+p(1));
end;


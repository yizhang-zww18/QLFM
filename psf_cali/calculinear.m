function [k,m] = calculinear(true_phase,calcu_phase);
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
% Y = kX+b;
Y = true_phase;
X = calcu_phase;
[rr,cc] = size(X);
ra = (rr-1)/2;
[xx,yy]=meshgrid([-ra:ra],[-ra:ra]);
mask = xx.^2+yy.^2<=(ra^2);
Yline = Y(ind2sub([rr,cc],find(mask)));
Xline = X(ind2sub([rr,cc],find(mask)));
SUMxy = sum(sum(Yline.*Xline));
SUMy = sum(Yline);
SUMx = sum(Xline);
SUMx2 = sum(Xline.^2);
SUM2x = sum(Xline).^2;
N = length(Xline);
k = (N*SUMxy-SUMy*SUMx)/(N*SUMx2-SUM2x);
m = (SUMy-k*SUMx)/N;
tt = 0;
end


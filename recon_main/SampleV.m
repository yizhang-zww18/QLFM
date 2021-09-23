function [UP_Xguess] = SampleV(Xguess,ind1,ind2)
% upsample 3D volume
% Xguess: the original volume
% ind1: original index along z axis
% ind2: upsampled index along z axis


[img_r,img_c,z] = size(Xguess);
x_index=1:size(Xguess,1);
y_index=1:size(Xguess,2);

[X1,Y1,Z1]=meshgrid(y_index,x_index,ind1);
[X2,Y2,Z2]=meshgrid(y_index,x_index,ind2);

X=gather(Xguess);
Y=interp3(X1,Y1,Z1,X,X2,Y2,Z2);
Y(Y<0) = 0;
Y(isnan(Y)) = 0;
UP_Xguess = Y;
end


function z=SH(c,data)
x=data(1,:);
y=data(2,:);
[X,Y]=meshgrid(x,y);
[theta,r]=cart2pol(X,Y);
idx=r<=1;
zz=zeros(size(X));
%mask = zeros(size(c));
%mask(7:end) = 1;
%c = c.*mask;
b = c(1:4);
a = c(5:end);
% a=c;

% zz(idx) = a*zernfun(4,0,r(idx),theta(idx));
% 
zz(idx)=b(1)*zernfun(0,0,r(idx),theta(idx))+b(2)*zernfun(1,1,r(idx),theta(idx))+b(3)*zernfun(1,-1,r(idx),theta(idx))+b(4)*zernfun(2,0,r(idx),theta(idx))+...
a(1)*zernfun(2,2,r(idx),theta(idx))+a(2)*zernfun(2,-2,r(idx),theta(idx))+a(3)*zernfun(3,1,r(idx),theta(idx))+...
a(4)*zernfun(3,-1,r(idx),theta(idx))+a(5)*zernfun(3,3,r(idx),theta(idx))+a(6)*zernfun(3,-3,r(idx),theta(idx))+...
a(7)*zernfun(4,0,r(idx),theta(idx))+a(8)*zernfun(4,2,r(idx),theta(idx))+a(9)*zernfun(4,-2,r(idx),theta(idx))+...
a(10)*zernfun(4,4,r(idx),theta(idx))+a(11)*zernfun(4,-4,r(idx),theta(idx))+a(12)*zernfun(5,1,r(idx),theta(idx))+...
a(13)*zernfun(5,-1,r(idx),theta(idx))+a(14)*zernfun(5,3,r(idx),theta(idx))+a(15)*zernfun(5,-3,r(idx),theta(idx))+...
a(16)*zernfun(5,5,r(idx),theta(idx))+a(17)*zernfun(5,-5,r(idx),theta(idx))+a(18)*zernfun(6,0,r(idx),theta(idx))+...
a(19)*zernfun(6,2,r(idx),theta(idx))+a(20)*zernfun(6,-2,r(idx),theta(idx))+a(21)*zernfun(6,4,r(idx),theta(idx))+...
a(22)*zernfun(6,-4,r(idx),theta(idx))+a(23)*zernfun(6,6,r(idx),theta(idx))+a(24)*zernfun(6,-6,r(idx),theta(idx))...
+a(25)*zernfun(7,1,r(idx),theta(idx))+a(26)*zernfun(7,-1,r(idx),theta(idx))+a(27)*zernfun(7,3,r(idx),theta(idx))...
+a(28)*zernfun(7,-3,r(idx),theta(idx))+a(29)*zernfun(7,5,r(idx),theta(idx))+a(30)*zernfun(7,-5,r(idx),theta(idx))...
+a(31)*zernfun(7,7,r(idx),theta(idx))+a(32)*zernfun(7,-7,r(idx),theta(idx))...
+a(33)*zernfun(8,0,r(idx),theta(idx))+a(34)*zernfun(8,2,r(idx),theta(idx))+a(35)*zernfun(8,-2,r(idx),theta(idx))...
+a(36)*zernfun(8,4,r(idx),theta(idx))+a(37)*zernfun(8,-4,r(idx),theta(idx))...
+a(38)*zernfun(8,6,r(idx),theta(idx))+a(39)*zernfun(8,-6,r(idx),theta(idx))...
+a(40)*zernfun(8,8,r(idx),theta(idx))+a(41)*zernfun(8,-8,r(idx),theta(idx));


z=zz;
% clc;
% clear;
load('../psf_shift/12102018_shift_cap180_with_sim_ite3.mat');

% %%%%% defocus %%%%%
% Nnum = 13;
% shift2 = shift_kernel;
% [Sx,Sy]=meshgrid([-fix(Nnum/2):fix(Nnum/2)],[-fix(Nnum/2):fix(Nnum/2)]);
% mask = (Sx.^2+Sy.^2)<=3.5^2;
% Sx=gather(Sx.*mask);
% Sy=gather(Sy.*mask);
% 
% 
% cx=shift2(round(Nnum/2),round(Nnum/2),1);
% cy=shift2(round(Nnum/2),round(Nnum/2),2);
% 
% shift = shift2;
% shift(:,:,1)=squeeze(shift_kernel(:,:,1))-cx;
% shift(:,:,2)=squeeze(shift_kernel(:,:,2))-cy;
% k1 = Sy.*squeeze(shift(:,:,1)).*mask+Sx.*squeeze(shift(:,:,2)).*mask;
% k2 = (Sx.*Sx+Sy.*Sy).*mask;
% k=sum(k1(:))/sum(k2(:));
% shift(:,:,1)=squeeze(shift(:,:,1))-k*Sy;
% shift(:,:,2)=squeeze(shift(:,:,2))-k*Sx;
% 
%%%%%%%%%%%%%%%%%%%%%%%
% shift = shift_kernel;
% shift_kernel(:,:,1) = shift(:,:,1)-shift(7,7,1);
% shift_kernel(:,:,2) = shift(:,:,2)-shift(7,7,2);
% shift_kernel(:,:,1) = shift_kernel(:,:,1) - shift_kernel(7,7,1);
% shift_kernel(:,:,2) = shift_kernel(:,:,2) - shift_kernel(7,7,2);

Nnum = 13;
[Sx,Sy]=meshgrid([-fix(Nnum/2):fix(Nnum/2)],[-fix(Nnum/2):fix(Nnum/2)]);
mask = (Sx.^2+Sy.^2)<=7^2;
shift_kernel(:,:,1) = shift_kernel(:,:,1).*mask;
shift_kernel(:,:,2) = shift_kernel(:,:,2).*mask;
% % 
% shift_kernel(:,:,1) = rot90(shift_kernel(:,:,1),2);
% shift_kernel(:,:,2) = rot90(shift_kernel(:,:,2),2);

waveShape = -shift_kernel;
% waveShape = -shift;
waveShape(abs(waveShape)>50) = 0;
[Nnum,~] = size(waveShape);
r_actual = 10.6;
expand = 5;
waveShape_expand = zeros(expand*Nnum,expand*Nnum,2);
for idu = 1:Nnum
    for idv = 1:Nnum
        waveShape_expand((idu-1)*expand+1:idu*expand,(idv-1)*expand+1:idv*expand,1) = waveShape(idu,idv,1);
        waveShape_expand((idu-1)*expand+1:idu*expand,(idv-1)*expand+1:idv*expand,2) = waveShape(idu,idv,2);
    end
end
[xx,yy] = meshgrid(-(expand*size(waveShape,1)-1)/2:(expand*size(waveShape,1)-1)/2,...
    -(expand*size(waveShape,1)-1)/2:(expand*size(waveShape,1)-1)/2);
mask = xx.^2+yy.^2<=((expand*r_actual/2).^2);
waveShape_expand = waveShape_expand.*mask;
waveShape_expand = waveShape_expand((end+1)/2-round(expand*r_actual/2):(end+1)/2+round(expand*r_actual/2),...
    (end+1)/2-round(expand*r_actual/2):(end+1)/2+round(expand*r_actual/2),:);

% load('../create_abe/aberphase_443_2107_qiucha18.mat');
ps = 193;
ns = 1015;
% true_phase = t;

[x1,y1] = meshgrid(1:size(waveShape_expand,1),1:size(waveShape_expand,2));
[x2,y2] = meshgrid(linspace(1,size(waveShape_expand,1),ps),linspace(1,size(waveShape_expand,1),ps));

calcu_dephase = zeros(ps,ps,2);
calcu_dephase(:,:,1)  = interp2(x1,y1,waveShape_expand(:,:,1),x2,y2,'nearest');
calcu_dephase(:,:,2)  = interp2(x1,y1,waveShape_expand(:,:,2),x2,y2,'nearest');

maxIte = 500;
calcu_phase = intercircle(calcu_dephase,maxIte);
% calcu_phase_inv = intercircle_inv(calcu_dephase,maxIte);
% [k,m] = calculinear(true_phase,calcu_phase);
% [k2,m2] = calculinear(true_phase,calcu_phase_inv);

[rr,cc] = size(calcu_phase);
ra = (rr-1)/2;
[xx,yy]=meshgrid([-ra:ra],[-ra:ra]);
mask = xx.^2+yy.^2<=(ra^2);
calcu_phase_k = 4.28*calcu_phase.*mask;
% calcu_phase_k2 = k2*calcu_phase_inv.*mask+m2;

x=linspace(-1,1,size(calcu_phase_k,1));
y=linspace(-1,1,size(calcu_phase_k,2));
xy=[x;y];
a1=lsqcurvefit('SH',zeros(1,45),xy,calcu_phase_k);
% a2= lsqcurvefit('SH',zeros(1,45),xy,calcu_phase_k2);


% figure;
% subplot(121);
% imshow(true_phase,[]);
% subplot(122);
% imshow(calcu_phase_k,[]);
figure;
imshow(calcu_phase_k,[]);
% 
% figure;
% subplot(121);
% imshow(exp(1i.*true_phase),[]);
% subplot(122);
% imshow(exp(1i.*calcu_phase_k),[]);
figure;
imshow(exp(1i.*(calcu_phase_k)),[]);


% t = calcu_phase_k;
% tt=padarray(calcu_phase_k,[(ns-ps)/2,(ns-ps)/2]);
% aber_phase=exp(1i.*tt);
% a1(4) = 0;
% a1(11) = a1(11)+20;
% save('12032329_zernike_para_ite2_10k_with_sim.mat','a1');
% % % 
% % 
% % a1(4) = 0;
save('zernike_para_ite3_rot180.mat','a1');
% % % % 
% a1(4) = 0;
% ps = 399;
% ns = 2107;
% a1(1:4) = 0;
% a1(16:end) = 0;
m = linspace(-1,1,ps);
n = linspace(-1,1,ps);
mn = [m;n];
phase1=SH(a1,mn);
tt=padarray(phase1,[(ns-ps)/2,(ns-ps)/2]);
aber_phase=exp(1i.*tt);
% save('12210_aberphase_399_2107_correctphase_ite3_k4.28_rot180.mat','aber_phase','phase1','tt');
% 
% save('aberphase_388_2029_calcu_qiucha.mat','aber_phase','tt','t');



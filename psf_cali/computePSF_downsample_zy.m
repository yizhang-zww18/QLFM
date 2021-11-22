% compute PSF of 5 dimensions
% pixel_x,pixel_y,object_x,object_y,object_z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parameter initialization %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
NA =        0.5;
MLPitch =   100*1e-6;% microlens pitch size
Nnum =      13;% virtual pixel number per microlens
fml =       2100*1e-6;% microlens focal length
lambda =    525*1e-9;% laser light wave length
OSR =       3;% calculate times in one pixel
n =         1;% refractive index
%M =         NA*2*fml/MLPitch;% magnification || formula(13)
M = 22.56;
%zmax =      (3.6+0.001)*1e-6;%
%zmin =      (-3.6+0.001)*1e-6;%
%zspacing =  0.3e-6;%

zmax = (500+0.001)*1e-6;
zmin = (-500+0.001)*1e-6;
zspacing = 1*1e-6;

eqtol = 1e-10;% threshold
tol = 0.005; % threshold

k = 2*pi*n/lambda; % k
k0 = 2*pi*1/lambda; % k
d = fml;   % optical distance between the microlens and the sensor
ftl = 165e-3;        % focal length of tube lens
fobj = ftl/M;  % focal length of objective lens
fnum_obj = M/(2*NA); % f-number of objective lens (imaging-side) || not used in this code
fnum_ml = fml/MLPitch; % f-number of micrl lens || not used in this code

if mod(Nnum,2)==0,
    error(['Nnum should be an odd number']);
end
pixelPitch = MLPitch/Nnum; %% pitch of virtual pixels


% define object space
% test PSF on line on the NATIVE image plane
x1objspace = [0];
x2objspace = [0];
x3objspace = [zmin:zspacing:zmax];
% x3objspace = -20*1e-6;

objspace = ones(length(x1objspace),length(x2objspace),length(x3objspace));% discrete object space

IMGsize=Nnum*43;
IMGsize_L = Nnum*43;
% psf=zeros(IMGsize,IMGsize,Nnum,Nnum,length(x3objspace));

p3max = max(abs(x3objspace));% the nearset z
x1testspace = (pixelPitch/OSR)* [0:1: Nnum*OSR*300];% test image space on x1 direction
x2testspace = [0];   %
[psfLine] = calcPSFFT(p3max, fobj, NA, x1testspace, pixelPitch/OSR, lambda, d, M, n);% line PSF on the native image plane
outArea = find(psfLine<0.02);%
if isempty(outArea),
    error('Estimated PSF size exceeds the limit');
end
IMGSIZE_REF = ceil(outArea(1)/(OSR*Nnum));%% effective microlens
disp(['Size of PSF ~= ' num2str(IMGSIZE_REF) ' [microlens pitch]' ]);

IMG_HALFWIDTH = max( Nnum*(IMGSIZE_REF + 1), 2*Nnum);
disp(['Size of IMAGE = ' num2str(IMG_HALFWIDTH*2*OSR+1) 'X' num2str(IMG_HALFWIDTH*2*OSR+1) '' ]); % pixel number
x1space = (pixelPitch/OSR)*[-IMG_HALFWIDTH*OSR:1:IMG_HALFWIDTH*OSR]; % x space
x2space = (pixelPitch/OSR)*[-IMG_HALFWIDTH*OSR:1:IMG_HALFWIDTH*OSR]; % y space
x1length = length(x1space);%x
x2length = length(x2space);%y

x1MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2]; % total x space per ML
x2MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2]; % total y space per ML
x1MLdist = length(x1MLspace);%
x2MLdist = length(x2MLspace);%

validpts = find(objspace>eqtol);% find non-zero points
numpts = length(validpts);%
[p1indALL p2indALL p3indALL] = ind2sub( size(objspace), validpts);% index to subcripts
p1ALL = x1objspace(p1indALL)';% effective obj points x location
p2ALL = x2objspace(p2indALL)';% effective obj points y location
p3ALL = x3objspace(p3indALL)';% effective obj points z location

MLARRAY = calcML(fml, k0, x1MLspace, x2MLspace, x1space, x2space); % micro array phase mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% projection from points on z axis %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LFpsfWAVE_STACK = zeros(x1length, x2length, numpts);%
%psfWAVE_STACK = zeros(x1length, x2length, numpts);%
disp(['Start Calculating PSF...']);

centerPT = ceil(length(x1space)/2);
halfWidth =  Nnum*(IMGSIZE_REF + 0 )*OSR;%
centerArea = (  max((centerPT - halfWidth),1)   :   min((centerPT + halfWidth),length(x1space))     );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Compute Light Field PSFs (light field) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];% corresponding object space x1
x2objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];% corresponding object space x2
XREF = ceil(length(x1objspace)/2);
YREF = ceil(length(x1objspace)/2);
CP = ( (centerPT-1)/OSR+1 - halfWidth/OSR :1: (centerPT-1)/OSR+1 + halfWidth/OSR  );%
H_z = zeros( length(CP), length(CP), length(x1objspace), length(x2objspace));%


disp(['Computing PSFs (1/3)']);
for eachpt=1:numpts      %
    time_s = tic;
    p1 = p1ALL(eachpt); % object point #eachpt x
    p2 = p2ALL(eachpt);
    p3 = p3ALL(eachpt);
    
    IMGSIZE_REF_IL = ceil(IMGSIZE_REF*( abs(p3)/p3max));%
    halfWidth_IL =  max(Nnum*(IMGSIZE_REF_IL + 0 )*OSR, 2*Nnum*OSR);
    centerArea_IL = (  max((centerPT - halfWidth_IL),1)   :   min((centerPT + halfWidth_IL),length(x1space))     );% corresponding center area
    disp(['point #',num2str(eachpt),' on the z axis:']);
    disp(['size of center area = ' num2str(length(centerArea_IL)) 'X' num2str(length(centerArea_IL)) ]);
    
    % excute PSF computing funcion
    [psfWAVE, LFpsfWAVE] = calcPSF(p1, p2, p3, fobj, NA, x1space, x2space, pixelPitch/OSR, lambda, MLARRAY, fml, M, n,  centerArea_IL);
    %psfWAVE_STACK(:,:,eachpt)  = psfWAVE;% PSF on the NATIVE image plane
    %LFpsfWAVE_STACK(:,:,eachpt)= LFpsfWAVE; %  PSF on image plane
    
    psfREF = psfWAVE;
    for b1 = 1:length(x2objspace)
        for a1 = 1:length(x1objspace)
            psfSHIFT0= im_shift2(psfREF, OSR*(a1-XREF), OSR*(b1-YREF) );%
            [f1,dx1,x1]=fresnel2D(psfSHIFT0.*MLARRAY, pixelPitch/OSR, fml,lambda);%
            f1= im_shift2(f1, -OSR*(a1-XREF), -OSR*(b1-YREF) );%
            f1_AP = f1;
            %f1_AP( (xmin:xmax), (ymin:ymax) ) = f1( (xmin:xmax), (ymin:ymax) );
            [f1_AP_resize, x1shift, x2shift] = pixelBinning(abs(f1_AP).^2, OSR);
            f1_CP = f1_AP_resize( CP - x1shift, CP-x2shift );
            H_z(:,:,a1,b1) = f1_CP;%
            ttt = 1;
        end
    end
    H4Dslice = H_z;
    H4Dslice(find(H4Dslice< (tol*max(H4Dslice(:))) )) = 0;% remove noise
    H_z = H4Dslice;
    
    disp('normalize');
    
    for b2 = 1:length(x2objspace)
        for a2 = 1:length(x1objspace)
            sss = H_z(:,:,a2,b2);
            H_z(:,:,a2,b2) = H_z(:,:,a2,b2)./sum(sss(:));
        end
    end
    
    disp('split');
    
    border=fix(IMGsize_L/2)-fix(size(H_z,1)/2);
    blur_image=zeros(IMGsize_L,IMGsize_L,size(H_z,3),size(H_z,4));
    for i=1:size(H_z,3)
        for j=1:size(H_z,4)
            temp=zeros(IMGsize_L,IMGsize_L);
            temp(border+1:end-border,border+1:end-border)=squeeze(H_z(:,:,i,j));
            blur_image(:,:,i,j)=im_shift3d(temp,i-((Nnum+1)/2),j-((Nnum+1)/2));
        end
    end
    
    bb=zeros(Nnum,Nnum,size(blur_image,1)/size(H_z,3),size(blur_image,2)/size(H_z,4),Nnum,Nnum);
    
    for i=1:size(H_z,3)
        for j=1:size(H_z,4)
            for a=1:size(blur_image,1)/size(H_z,3)
                for b=1:size(blur_image,2)/size(H_z,4)
                    bb(i,j,a,b,:,:)=squeeze(  blur_image(  (a-1)*Nnum+i,(b-1)*Nnum+j,:,:  )  );
                end
            end
        end
    end
    
    bbb=zeros(  size(blur_image,1),size(blur_image,2),Nnum,Nnum  );
    for a=1:size(blur_image,1)/size(H_z,3)
        for c=1:Nnum
            x=Nnum*a+1-c;
            for b=1:size(blur_image,2)/size(H_z,4)
                for d=1:Nnum
                    y=Nnum*b+1-d;
                    bbb(x,y,:,:)=squeeze(bb(:,:,a,b,c,d));
                end
            end
        end
    end
    
    %     for idu = 1:Nnum
    %         for idv = 1:Nnum
    %             nr = round(IMGsize_L/downsample_xy(eachpt));
    %             if(mod(nr,2)==0)
    %                 nr = nr+1;
    %             end
    %             ttt = squeeze(bbb(:,:,idu,idv));
    %             sss = sum(ttt(:));
    %             tmp = imresize(ttt,[nr,nr],'nearest');
    %             tmp2 = zeros(IMGsize_L,IMGsize_L);
    %             tmp2((IMGsize_L+1)/2-(nr-1)/2:(IMGsize_L+1)/2+(nr-1)/2,(IMGsize_L+1)/2-(nr-1)/2:(IMGsize_L+1)/2+(nr-1)/2) = tmp;
    %             if(sum(tmp2(:))>0)
    %                 tmp2 = tmp2./sum(tmp2(:))*sss;
    %             end
    %             bbb(:,:,idu,idv)=tmp2;
    %         end
    %     end
    
    tmp = bbb((IMGsize_L+1)/2-(IMGsize-1)/2:(IMGsize_L+1)/2+(IMGsize-1)/2,(IMGsize_L+1)/2-(IMGsize-1)/2:(IMGsize_L+1)/2+(IMGsize-1)/2,:,:);
    psf(:,:,:,:,eachpt) = tmp;
    disp(['idz = ',num2str(eachpt)]);
    time_e = toc(time_s);
    disp(['calcu z = ',num2str(eachpt),' layer take time ',num2str(time_e),' sec......']);
    temp=gather(psf(:,:,:,:,1));
    Output=uint16(gather(temp)./0.005.*65535);
    for idu=1:13%�����ؽ�ͼ��
        for idv = 1:13
            ii = (idu-1)*13+idv;
            if ii==1
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat('psf_sim','.tif'));
            else
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat('psf_sim','.tif'),'WriteMode', 'append');
            end
        end
    end
end
% psf2 = single(psf);
% save('psf_23x_inter_5_10.mat','psf2','-v7.3');


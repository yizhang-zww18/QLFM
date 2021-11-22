clc;
clear;
gpuDevice(2);
load('H:\project\02_edofLF\code\20210415_simu_10x_0.25NA_psf\genephase\20210420_aber_phase_nodefocus_513_2431.mat');
NA =        0.28;
MLPitch =   100*1e-6;   % microlens pitch size
Nnum =      15;         % virtual pixel number per microlens
fml =       1800*1e-6;  % microlens focal length
lambda =    590*1e-9;   % laser light wave length
OSR =       3;          % calculate times in one pixel
n =         1;          % refractive index
%M =         NA*2*fml/MLPitch;% magnification || formula(15)
M = 10;

zmax = (10+0.001)*1e-6;
zmin = (10+0.001)*1e-6;
zspacing = 5e-6;

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
% x3objspace = -260*1e-6;

%
% drate = [1,2,4,8,16,20];
% dlayer = [101,100,200,400,800,200];
% ind = [];
% block = length(drate);
% lengthZ = length(x3objspace);
% clp = (lengthZ+1)/2;
% for i = block:-1:1
%     if(i==block)
%         ind = 1:drate(i):dlayer(i)/2-1;
%     else
%         ee = ind(end);
%         ind = [ind,ee+drate(i+1):drate(i):ee+drate(i+1)+floor(dlayer(i)/2)-1];
%     end
% end
% ind = [ind,clp];
% for i = 1:block
%     ee = ind(end);
%     if(i==1)
%         ind = [ind,ee+1:ee+floor(dlayer(i)/2)];
%     else
%         ind = [ind,ee+drate(i):drate(i):ee+drate(i)+dlayer(i)/2-1];
%     end
% end
%
% save('psf_select_index.mat','ind');
% ind_z = ind;
% ind = ind-clp;
% ind = ind*zspacing;
% ind = ind+0.001*1e-6;
% x3objspace = ind;
objspace = ones(length(x1objspace),length(x2objspace),length(x3objspace));% discrete object space

IMGsize=Nnum*55;
IMGsize_L = Nnum*55;
HALF_ML_NUM = 26;
% psf=zeros(IMGsize,IMGsize,Nnum,Nnum,length(x3objspace));

validpts = find(objspace>eqtol);% find non-zero points
numpts = length(validpts);%
[p1indALL p2indALL p3indALL] = ind2sub( size(objspace), validpts);% index to subcripts
p1ALL = x1objspace(p1indALL)';% effective obj points x location
p2ALL = x2objspace(p2indALL)';% effective obj points y location
p3ALL = x3objspace(p3indALL)';% effective obj points z location


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% projection from points on z axis %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LFpsfWAVE_STACK = zeros(x1length, x2length, numpts);%
%psfWAVE_STACK = zeros(x1length, x2length, numpts);%
disp(['Start Calculating PSF...']);

% centerPT = ceil(length(x1space)/2);
% halfWidth =  Nnum*(IMGSIZE_REF + 0 )*OSR;%
% centerArea = (  max((centerPT - halfWidth),1)   :   min((centerPT + halfWidth),length(x1space))     );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Compute Light Field PSFs (light field) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%% ∂‘±» defocus phase ∫Õ ÷±Ω”À„ %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixelPitch_OSR = MLPitch/OSR/Nnum; %simulated pixel size after OSR
fov=length(-(HALF_ML_NUM+1)*OSR*Nnum:(HALF_ML_NUM+1)*OSR*Nnum)*pixelPitch_OSR;   %the size of field of view for the PSF
pixelSize_OSR=length(-(HALF_ML_NUM+1)*OSR*Nnum:(HALF_ML_NUM+1)*OSR*Nnum); %the number of the pixels for each PSF
k2=2*pi/lambda;

sinalpha_max = NA / n / M;
% fx_sinalpha = lambda / 2 / pixelPitch2 / lambda;
fx_sinalpha = 1/(2*pixelPitch_OSR);
fx_step = 1/fov ;
fx_max = fx_sinalpha ;
fx= -fx_max+fx_step/2 : fx_step : fx_max;
[fxcoor fycoor] = meshgrid( fx , fx );
fx2coor=fxcoor.*fxcoor;
fy2coor=fycoor.*fycoor;

aperture_mask=((fx2coor+fy2coor)<=((NA/(lambda*M)).^2));
psfWAVE2=ones(pixelSize_OSR,pixelSize_OSR).*aperture_mask;
psfWAVE2 = gpuArray(single(psfWAVE2));
% tmp = fftshift(fft2(psfWAVE2));


x1MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2]; % total x space per ML
x2MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2]; % total y space per ML
x1MLdist = length(x1MLspace);%
x2MLdist = length(x2MLspace);%

x1space = (pixelPitch/OSR)*[-(HALF_ML_NUM+1)*Nnum*OSR:1:(HALF_ML_NUM+1)*Nnum*OSR]; % x space
x2space = (pixelPitch/OSR)*[-(HALF_ML_NUM+1)*Nnum*OSR:1:(HALF_ML_NUM+1)*Nnum*OSR]; % y space
x1length = length(x1space);%x
x2length = length(x2space);%y

MLARRAY = calcML(fml, k0, x1MLspace, x2MLspace, x1space, x2space); % micro array phase mask
MLARRAY = gpuArray(single(MLARRAY));

x1objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];% corresponding object space x1
x2objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];% corresponding object space x2
XREF = ceil(length(x1objspace)/2);
YREF = ceil(length(x1objspace)/2);
centerPT = ceil(length(x1space)/2);
halfWidth = HALF_ML_NUM*Nnum*OSR;
CP = ( (centerPT-1)/OSR+1 - halfWidth/OSR :1: (centerPT-1)/OSR+1 + halfWidth/OSR  );%
H_z = zeros(length(CP),length(CP),15,15);
% H = zeros(length(CP),length(CP),15,15,length(x3objspace));
psf_z = zeros(IMGsize,IMGsize,15,15,'single');
% psf_z = zeros(IMGsize,IMGsize,15,15,'single');
% H_z = gpuArray(single(H_z));
% load('../create_abe/aberphase_147_781.mat','aber_phase');
filepath = strcat('../10x_0.28na_aph3neg6.7_neg300T300_dz5');
% filepath = strcat('../test');
mkdir(filepath);
% aber = ones(size(psfWAVE2));
% for iteshift = [0:2]
%     load(strcat('../zernike/aberphase_147_781_cap_sim_calcu_k0.03_shift',num2str(iteshift),'.mat'),'aber_phase');
%     aber = aber.*aber_phase;
% end
% aber_phase = aber;
% load('../create_abe/aberphase_185_859_qiucha5.mat','aber_phase');
% aber1 = aber_phase;
% load('../zernike/12210_aberphase_399_2107_correctphase_ite3_k4.28_rot180.mat','aber_phase');
% aber2 = aber_phase;
% aber2 = rot90(aber2,3);
aber_phase = aber_phaseV_exp;
for eachpt=1:numpts
    aa = tic;
    if(eachpt<0)
        continue;
    else
        disp(['calcu point #',num2str(eachpt),' ...............']);
        time_s = tic;
        p1 = p1ALL(eachpt); % object point #eachpt x
        p2 = p2ALL(eachpt);
        p3 = p3ALL(eachpt);
        
        timeWAVE = tic;
        tempP=k2*n*p3*realsqrt((1-(fxcoor.*lambda./n.*M).^2-(fycoor.*lambda./n.*M).^2).*aperture_mask);
        tempP = gpuArray(single(tempP));
        psfWAVE_fAFTERNO=psfWAVE2.*exp(1j*tempP);
        psfWAVE_AFTERNO=fftshift(ifft2(ifftshift(squeeze(psfWAVE_fAFTERNO))));
        psfWAVE_AFTERNO = ifft2(ifftshift(fftshift(fft2(psfWAVE_AFTERNO)).*aber_phase));
        %         psfWAVE_AFTERNO = ifft2(ifftshift(fftshift(fft2(psfWAVE_AFTERNO)).*aber2));
        tt = toc(timeWAVE);
        disp(['calcu psfWAVE take ',num2str(tt),' sec....']);
        
        %     psfWAVE_AFTERNO=abs(squeeze(psfWAVE_AFTERNO));
        
        timeFre = tic;
        %     psfREF = gpuArray(single(psfWAVE_AFTERNO));
        %     psfREF = psfWAVE_AFTERNO;
        for b1 = 1:length(x2objspace)
            for a1 = 1:length(x1objspace)
                timein = tic;
                psfSHIFT0= im_shift2_GPU(psfWAVE_AFTERNO, OSR*(a1-XREF), OSR*(b1-YREF) );%
                f1=fresnel2D_GPU(psfSHIFT0.*MLARRAY, pixelPitch/OSR, fml,lambda);%
                f1= im_shift2_GPU(f1, -OSR*(a1-XREF), -OSR*(b1-YREF) );%
                [f1_AP_resize, x1shift, x2shift] = pixelBinning_GPU(abs(f1).^2, OSR);
                f1_CP = f1_AP_resize( CP - x1shift, CP-x2shift );
                H_z(:,:,a1,b1) = gather(f1_CP);%
                tt = toc(timein);
                disp(['calcu one point H take ',num2str(tt),' sec....']);
            end
        end
        tt = toc(timeFre);
        disp(['calcu H take ',num2str(tt),' sec....']);
        
        H4Dslice = H_z;
        H4Dslice(find(H4Dslice< (tol*max(H4Dslice(:))) )) = 0;% remove noise
        H_z = H4Dslice;
        %         H(:,:,:,:,eachpt) = H_z;
        
        disp(['normalize...NA_',num2str(NA)]);
        
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
                blur_image(:,:,i,j)=(im_shift3d(temp,i-((Nnum+1)/2),j-((Nnum+1)/2)));
                %                 blur_image(:,:,i,j) = temp;
            end
        end
        
        blur_image(isnan(blur_image)) = 0;
        maxH_z  = max(blur_image(:));
        Output=uint16(blur_image./maxH_z*65535);
        for idu=1:15
            for idv = 1:15
                ii = (idu-1)*15+idv;%ÔøΩÔøΩÔøΩÔøΩÔøΩÿΩÔøΩÕºÔøΩÔøΩ
                if ii==1
                    imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat('H_z','.tif'));
                else
                    imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat('H_z','.tif'),'WriteMode', 'append');
                end
            end
        end
        command = sprintf('ReAlign %d %s %s %d %d %d %d %d %d %s %d %d %d %d',...
            15, 'psf_temp', 'H_z.tif', 0, (IMGsize_L-1)/2, (IMGsize_L-1)/2, ((IMGsize_L/15)-1)/2,...
            1, 0, '15x15.conf.sk.png',1,0, 15, ((IMGsize_L/15)-1)/2);
        system(command);
        
        bbb = zeros(IMGsize_L,IMGsize_L,15,15);
        for idu = 1:15
            for idv = 1:15
                tmp = imread('psf_temp_No0.tif',(idu-1)*15+idv);
                bbb(:,:,idu,idv) = double(tmp)./65535.*maxH_z;
            end
        end
        
        tmp = bbb((IMGsize_L+1)/2-(IMGsize-1)/2:(IMGsize_L+1)/2+(IMGsize-1)/2,(IMGsize_L+1)/2-(IMGsize-1)/2:(IMGsize_L+1)/2+(IMGsize-1)/2,:,:);
        %         tmp = flip(tmp,1);
        %         tmp = flip(tmp,2);
        psf_z = single(tmp);
        temp=gather(tmp);
        Output=uint16(gather(temp)./0.01.*65535);
        for idu=1:15%ÔøΩÔøΩÔøΩÔøΩÔøΩÿΩÔøΩÕºÔøΩÔøΩ
            for idv = 1:15
                ii = (idu-1)*15+idv;
                if ii==1
                    imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(filepath,'/psf_sim',num2str(eachpt),'.tif'));
                else
                    imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(filepath,'/psf_sim',num2str(eachpt),'.tif'),'WriteMode', 'append');
                end
            end
        end
        
        save([filepath,'/psf_sim_10_',num2str(eachpt),'.mat'],'psf_z','-v7.3');
        onetime = toc(aa);
        disp(['idz = ',num2str(eachpt),'taketime ',num2str(onetime),' sec......']);
        
    end
end
% save(strcat(filepath,'/psf10_qiucha10.5_neg15_15',num2str(NA),'.mat'),'psf','-v7.3');
% save('../deAber_sim/psf_neg5_0_pos5.mat','psf');
% filepath = './save_psf_23x_z50_dz1';
% mkdir(filepath);
% save(strcat(filepath,'/psf_23x_z50_dz1.mat'),'psf','-v7.3');
% time_e = toc(time_s);
% disp(['calcu z = ',num2str(eachpt),' layer take time ',num2str(time_e),' sec......']);

% psf2 = single(psf);
% save('psf_23x_defocus_5_10.mat','psf2','-v7.3');


clc;
clear;
iteration = 20   ;
gpuDevice(1);
M =         5.6509; % M = NA*2*fml/MLPitch;% magnification || formula(15) ?
Nnum =      15; % 15 virtual pixels per microlens
n =         1; % refractive index
NA =        0.5; % NA of objective lens is 0.5.
MLPitch =   69*1e-6; % The pitch size of microlens is 69 um
fml =       393.3*1e-6; % microlens focal length 193 um
ftl =       265e-3; % Focal length of tube lens is 165 mm
OSR =       3; % calculate 5 times per pixel
lambda =    525*1e-9; %  light wave length is 525 nm

zmax = (75+0.001)*1e-6;
zmin = (75+0.001)*1e-6;
zspacing = 1e-6;

eqtol = 1e-10;% threshold
tol = 0.005; % threshold

k =         2*pi*n/lambda;  % k
k0 =        2*pi*1/lambda;  % k
d =         fml;            % optical distance between the microlens and the sensor
% ftl =       200e-3;         % focal length of tube lens
fobj =      ftl/M;          % focal length of objective lens
fnum_obj =  M/(2*NA);       % f-number of objective lens (imaging-side) || not used in this code
fnum_ml =   fml/MLPitch;    % f-number of micrl lens || not used in this code

if mod(Nnum,2)==0,
    error(['Nnum should be an odd number']);
end
pixelPitch = MLPitch/Nnum; %% pitch of virtual pixels

% define object space
% test PSF on line on the NATIVE image plane

x1objspace = [0];
x2objspace = [0];
x3objspace = [zmin:zspacing:zmax];
objspace = ones(length(x1objspace),length(x2objspace),length(x3objspace));% discrete object space

IMGsize=Nnum*55;
IMGsize_L = Nnum*55;
HALF_ML_NUM = 26;


date = '2021_1122_air_pos75_rot_onlysphe';
phaseangle = ['../',date,'_phasePic'];
mkdir(phaseangle);
filepath = strcat('../',date,'_6x_psf_sim');
mkdir(filepath);
shiftpath = strcat('../',date,'_psf_shift');
mkdir(shiftpath);
totalPhase = zeros(1257,1257,11);
totalInten = ones(2431,2431,11);

for ite = 1:iteration
    %% calcu psf
    validpts = find(objspace>eqtol);% find non-zero points
    numpts = length(validpts);%
    [p1indALL p2indALL p3indALL] = ind2sub( size(objspace), validpts);% index to subcripts
    p1ALL = x1objspace(p1indALL)';% effective obj points x location
    p2ALL = x2objspace(p2indALL)';% effective obj points y location
    p3ALL = x3objspace(p3indALL)';% effective obj points z location
    disp(['Start Calculating PSF...']);
    
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
    
    % aber = ones(size(psfWAVE2));
    % for iteshift = [0:2]
    %     load(strcat('../zernike/aberphase_147_781_cap_sim_calcu_k0.03_shift',num2str(iteshift),'.mat'),'aber_phase');
    %     aber = aber.*aber_phase;
    % end
    % aber_phase = aber;
    % load('../create_abe/aberphase_185_859_qiucha5.mat','aber_phase');
    % aber1 = aber_phase;
    
    
    for eachpt=1:numpts
        aa = tic;
        if(eachpt<0)
            continue;
        else
%         elseif(eachpt==17)
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
            if(ite~=1)
                load([phaseangle,'/aberphase_1257_2431_layer_',num2str(eachpt),'_ite',num2str(ite-1),'.mat']);
%                 load([phaseangle,'/intensitymask_layer_',num2str(sim_layer),'_ite',num2str(ite-1),'.mat'],'intensity_mask');
                psfWAVE_AFTERNO = ifft2(ifftshift(fftshift(fft2(psfWAVE_AFTERNO)).*aber_phase));
            else
%                 load('../create_abe/aberphase_1257_2431_1iucha_neg3.mat');
%                 psfWAVE_AFTERNO = ifft2(ifftshift(fftshift(fft2(psfWAVE_AFTERNO)).*aber_phase));
%                 load(['../2019_1216_1401_phaseonly_phasePic/aberphase_1257_2431_layer_',num2str(eachpt),'_ite20.mat']);
%                 psfWAVE_AFTERNO = ifft2(ifftshift(fftshift(fft2(psfWAVE_AFTERNO)).*aber_phase));
            end
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
%                     f1 = im_shift2_GPU(f1,1,-1);
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
                    ii = (idu-1)*15+idv;%�����ؽ�ͼ��
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
            Output=uint16(gather(temp)./0.1.*65535);
            for idu=1:15%�����ؽ�ͼ��
                for idv = 1:15
                    ii = (idu-1)*15+idv;
                    if ii==1
                        imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(filepath,'/psf_sim',num2str(eachpt),'_ite',num2str(ite-1),'.tif'));
                    else
                        imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(filepath,'/psf_sim',num2str(eachpt),'_ite',num2str(ite-1),'.tif'),'WriteMode', 'append');
                    end
                end
            end
            
            save([filepath,'/psf_sim_10_',num2str(eachpt),'_ite',num2str(ite-1),'.mat'],'psf_z','-v7.3');
            onetime = toc(aa);
            disp(['idz = ',num2str(eachpt),'taketime ',num2str(onetime),' sec......']);
        end
    end
    
    %% calcu shift
    for sim_layer = 1
        load([filepath,'/psf_sim_10_',num2str(sim_layer),'_ite',num2str(ite-1),'.mat']);
        psf1 = psf_z;
        
        file = ['D:\lab\project\spherical_aberation\data\2021-11-22_psf\movecenterpos0.5_nosphe_z75_5x5_170.0ms_Full_Hardware_LaserCount1_211119210341.0.tiff_No0-1.tif'];
        Info=imfinfo(file);
        Slice=size(Info,1);
        Width=Info.Width;
        Height=Info.Height;
        psf2 = zeros(Width,Height,15,15);
        for idu = 1:15
            for idv = 1:15
                ii = (idu-1)*15+idv;
                tmp = imread(file,ii);
                psf2(:,:,idu,idv) = rot90(tmp,2);%% 可能换成下一行
%                 psf2(:,:,idu,idv) = tmp;
            end
        end
        %
        [r1,c1,~] = size(psf1);
        [r2,c2,~] = size(psf2);
        
        ra = min(r1,r2);
        if(mod(ra,2)==0)
            ra = ra-1;
        end
        psf1 = psf1(round((r1+1)/2)-(ra-1)/2:round((r1+1)/2)+(ra-1)/2,round((r1+1)/2)-(ra-1)/2:round((r1+1)/2)+(ra-1)/2,...
            :,:);
        psf2 = psf2(round((r2+1)/2)-(ra-1)/2:round((r2+1)/2)+(ra-1)/2,round((r2+1)/2)-(ra-1)/2:round((r2+1)/2)+(ra-1)/2,...
            :,:);
        psf2 = psf2-100;
        psf2(psf2<0) = 0;
        
        Nnum = 15;
        shift_kernel = zeros(Nnum,Nnum,2);
        [xx,yy] = meshgrid([-(ra-1):(ra-1)],[-(ra-1):(ra-1)]);
        mask = xx.^2+yy.^2<=(ra/4)^2;
        for idu = 1:Nnum
            for idv = 1:Nnum
                img1 = squeeze(psf1(:,:,idu,idv));
                img2 = squeeze(psf2(:,:,idu,idv));
                if(sum(img1(:))>0 && sum(img2(:))>0)
                    corr = normxcorr2(img1,img2);
                    corr = corr.*mask;
                    [is1,is2] = find(corr==max(corr(:)),1);
                    a = (1:size(corr,1))';
                    b = corr(is1,:)';
                    f1 = fit((1:size(corr,1))',corr(is1,:)','gauss1','Lower',[0,is2-20,-inf],'Upper',[1,is2+20,inf]);
                    f2 = fit((1:size(corr,1))',corr(:,is2),'gauss1','Lower',[0,is1-20,-inf],'Upper',[1,is1+20,inf]);
                    testa=f2.b1;
                    testb=f1.b1;
                    shift_kernel(idu,idv,1) = testa-ra;
                    shift_kernel(idu,idv,2) = testb-ra;
                end
            end
            disp(['idu = ',num2str(idu),' ....']);
        end
%         map1 = squeeze(sum(sum(psf1,1),2));
% %         map1 = squeeze(max(max(psf1,[],1),[],2));
% %         map1 = map1./max(map1(:));
%         map2 = squeeze(sum(sum(psf2,1),2));
% %         map2 = squeeze(max(max(psf2,[],1),[],2))-100;
% %         map2(map2<0) = 0;
%         map2 = map2./max(map2(:));
        save([shiftpath,'/shift_cap180_with_sim_layer_',num2str(sim_layer),'_ite',num2str(ite),'.mat'],'shift_kernel');
    end
    %% calcu phase
    for sim_layer = 1
        load([shiftpath,'/shift_cap180_with_sim_layer_',num2str(sim_layer),'_ite',num2str(ite),'.mat'],'shift_kernel');
        Nnum = 15;
        [Sx,Sy]=meshgrid([-fix(Nnum/2):fix(Nnum/2)],[-fix(Nnum/2):fix(Nnum/2)]);
        mask = (Sx.^2+Sy.^2)<=8^2;
        shift_kernel(:,:,1) = shift_kernel(:,:,1).*mask;
        shift_kernel(:,:,2) = shift_kernel(:,:,2).*mask;
        
        waveShape = -shift_kernel;
        % waveShape = -shift;
        waveShape(abs(waveShape)>50) = 0;
        [Nnum,~] = size(waveShape);
        r_actual = 14;
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
        
        ps = 101;
%         ns = 1015;
        
        [x1,y1] = meshgrid(1:size(waveShape_expand,1),1:size(waveShape_expand,2));
        [x2,y2] = meshgrid(linspace(1,size(waveShape_expand,1),ps),linspace(1,size(waveShape_expand,1),ps));
        
        calcu_dephase = zeros(ps,ps,2);
        calcu_dephase(:,:,1)  = interp2(x1,y1,waveShape_expand(:,:,1),x2,y2,'nearest');
        calcu_dephase(:,:,2)  = interp2(x1,y1,waveShape_expand(:,:,2),x2,y2,'nearest');
        
        maxIte = 1000;
        calcu_phase = intercircle_zy_v2(calcu_dephase,maxIte);
        
        [rr,cc] = size(calcu_phase);
        ra = (rr-1)/2;
        [xx,yy]=meshgrid([-ra:ra],[-ra:ra]);
        mask = xx.^2+yy.^2<=(ra^2);
        calcu_phase_k = 4.28*calcu_phase.*mask;
        
        x=linspace(-1,1,size(calcu_phase_k,1));
        y=linspace(-1,1,size(calcu_phase_k,2));
        xy=[x;y];
        a1=lsqcurvefit('SH',zeros(1,45),xy,calcu_phase_k);
        %
        %     figure;
        %     imshow(calcu_phase_k,[]);
        %
        %     figure;
        %     imshow(exp(1i.*(calcu_phase_k)),[]);
        
        save([phaseangle,'/zernike_para_layer_',num2str(sim_layer),'_ite',num2str(ite),'_rot180.mat'],'a1');
        % % % %
        
        ps = 1257;
        ns = 2431;
        mm = linspace(-1,1,ps);
        nn = linspace(-1,1,ps);
        mn = [mm;nn];
        
        phase1=SH(a1,mn);
        tmp = angle(exp(1i.* phase1));
        tmp = (tmp-min(tmp(:)))./(max(tmp(:))-min(tmp(:)));
        imwrite(tmp,strcat(phaseangle,'/phaseadd_layer_',num2str(sim_layer),'_ite',num2str(ite),'.png'),'png');
        
        a1(1:4) = 0;
        a1(16:end) = 0;
%         a1(16:end) = 0;
%         if(ite==1)
%             a1(11) = a1(11)-3;
%         end
%         a1(abs(a1)<1) = 0;
        phase1=SH(a1,mn);
        totalPhase(:,:,sim_layer) = totalPhase(:,:,sim_layer)+phase1;
        
        [x2,y2] = meshgrid(-(ps-1)/2:(ps-1)/2,-(ps-1)/2:(ps-1)/2);
        mask = x2.^2+y2.^2<=((ps-1)/2)^2;
        totalPhase(:,:,sim_layer) = totalPhase(:,:,sim_layer).*mask;
        
        tmp = angle(exp(1i.* totalPhase(:,:,sim_layer)));
        tmp = (tmp-min(tmp(:)))./(max(tmp(:))-min(tmp(:)));
        imwrite(tmp,strcat(phaseangle,'/phasetotal_layer_',num2str(sim_layer),'_ite',num2str(ite),'.png'),'png');
        
        tmp = totalPhase(:,:,sim_layer);
        tt=padarray(tmp,[(ns-ps)/2,(ns-ps)/2]);
        aber_phase=exp(1i.*tt);
        save([phaseangle,'/aberphase_1257_2431_layer_',num2str(sim_layer),'_ite',num2str(ite),'.mat'],'aber_phase');
%         
        
%         map = map2./map1;
%         [xs,ys] = meshgrid(-6:6,-6:6);
%         mask_small = xs.^2+ys.^2<=((r_actual-1)/2)^2;
%         map = map.*mask_small;
%         map = map+~mask_small;
%         
%         map(isnan(map)) = 0;
%         map(isinf(map)) = 0;
%         r_new = round(ps/r_actual*15);
%         if(mod(r_new,2)==0)
%             r_new = r_new+1;
%         end
% 
%         map_large = imresize(map,[r_new,r_new],'cubic');
%         map_large(isnan(map_large)) = 0;
%         map_large(isinf(map_large)) = 0;
%         map_large(map_large<0) = 0;
%         map_large = map_large((r_new+1)/2-(ps-1)/2:(r_new+1)/2+(ps-1)/2,(r_new+1)/2-(ps-1)/2:(r_new+1)/2+(ps-1)/2);
%         map_large = map_large.*mask;
%         
%         map_large = map_large./sum(map_large(:)).*sum(mask(:));
%         b1=lsqcurvefit('SH',zeros(1,45),mn,double(map_large));
%         b1([2:3,5:10,12:21,23:36,38:45])=0;
%         map_large = SH(b1,mn);
%         map_large = map_large.*mask;
%         map_large = map_large./sum(map_large(:)).*sum(mask(:));
%         map_large = map_large+~mask;
% %         map_large = map_large.*mask;
%         map_large(isnan(map_large)) = 0;
%         map_large(isinf(map_large)) = 0;
%         map_large(map_large<0) = 0;
%         map_large = padarray(map_large,[(ns-ps)/2,(ns-ps)/2],1);
% %         map_large(isinf(map_large)) = 0;
%         tmp = map_large;
%         tmp = (tmp-min(tmp(:)))./(max(tmp(:))-min(tmp(:)));
%         imwrite(tmp,strcat(phaseangle,'/intensity_layer_',num2str(sim_layer),'_ite',num2str(ite),'.png'),'png');
%         save([phaseangle,'/aberphase_1257_2431_layer_',num2str(sim_layer),'_ite',num2str(ite),'.mat'],'aber_phase');
%         
%         [x3,y3] = meshgrid(-(ns-1)/2:(ns-1)/2,-(ns-1)/2:(ns-1)/2);
%         mask_large = x3.^2+y3.^2<=((ps-1)/2)^2;
%         stmp = totalInten(:,:,sim_layer).*map_large.*mask_large;
% %         stmp = tmp.*mask_large;
%         totalInten(:,:,sim_layer) = stmp./(sum(stmp(:)))*sum(mask_large(:))+~mask_large;
%         
%         tmp = totalInten(:,:,sim_layer);
%         intensity_mask = tmp;
%         tmp = (tmp-min(tmp(:)))./(max(tmp(:))-min(tmp(:)));
%         imwrite(tmp,strcat(phaseangle,'/intensity_total_layer_',num2str(sim_layer),'_ite',num2str(ite),'.png'),'png');
%         save([phaseangle,'/intensitymask_layer_',num2str(sim_layer),'_ite',num2str(ite),'.mat'],'intensity_mask');
        
        disp(['ite',num2str(ite),'_calcu_finished...........']);
    end
end
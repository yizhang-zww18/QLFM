%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processes and 3D reconstruciton with quantitative light field microscopy
% The Code is created based on the method described in the following paper
%   [1] Computational optical sectioning with an incoherent multiscale scattering model for light-field microscopy
%
%    Contact: YI ZHANG (zyi18@mails.tsinghua.edu.cn)
%    Date  : 09/23/2021

clc;
clear;
Methods = 'QLFM_sca_2D';% DECONVOLUTION MODEL.optional parameters: 'QLFM_ideal','QLFM_cali','QLFM_sca_LTR','QLFM_sca_RTL','QLFM_sca_2D';

switch Methods
    case 'QLFM_ideal'
        psf_file = '../psf/psf_ideal2.mat';
        savepath = 'Recon_QLFMwithoutScatter_ideal';
    case 'QLFM_cali'
        psf_file = '../psf/psf_cali2.mat';
        savepath = 'Recon_QLFMwithoutScatter_cali';
    case 'QLFM_sca_LTR'
        psf_file = '../psf/psf_cali2.mat';
        savepath = 'Recon_QLFMwithScatter_Left2Right';
    case 'QLFM_sca_RTL'
        psf_file = '../psf/psf_cali2.mat';
        savepath = 'Recon_QLFMwithScatter_Right2Left';
    case 'QLFM_sca_2D'
        psf_file = '../psf/psf_cali2.mat';
        savepath = 'Recon_QLFMwithScatter_2directions';
    otherwise
        msg = 'No method found';
        error(msg);
end

gpuDevice(1);% GPU accelerator
Nnum=13;% the number of sub-aperture components in one direction
load(psf_file);
load('../psf/index.mat');
[psf_r,psf_c,allu,allz] = size(psf);
index_2 = index_1(1):2:index_1(end);% multiscale index
mkdir(savepath);
file = '../dataset/c1-1.tif'; % load data
Info=imfinfo(file);
Slice=size(Info,1);
Width=Info.Width;
Height=Info.Height;
img_r = Height;
img_c = Width;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% reconstruction %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

func = @(x)8.916e-09.*x.^4-3.246e-06.*x.^3+ 0.0004943.*x.^2 -0.0362.*x+1.036; %zupdate
maxIter =2; % maximum number of iterations
UPweight = 0.6;% update rate
alpha = 1e-9; % learning rate
Vthres = 1;% threshold for V
CUTXY_BORDER = 50;
CUTZ_BORDER = 150;


if(psf_r>img_r)
    psf = psf((psf_r+1)/2-(img_r-1)/2:(psf_r+1)/2+(img_r-1)/2,:,:,:);
end
if(psf_c>img_c)
    psf = psf(:,(psf_c+1)/2-(img_c-1)/2:(psf_c+1)/2+(img_c-1)/2,:,:);
end

for frame = 0
    blur_image = zeros(img_r,img_c,13,13);
    for idu = 1:13
        for idv = 1:13
            tmp = imread(file,((idu-1)*13+idv));
            blur_image(:,:,idu,idv) = rot90(tmp,2);
        end
    end
    blur_image = blur_image-100;
    blur_image(blur_image<=0) = 0;
    
    Xguess=ones(size(blur_image,1),size(blur_image,2),allz);
    zupdate = zeros(size(Xguess));
    for idz = 1:181
        zupdate(:,:,idz) = ones(size(blur_image,1),size(blur_image,2)).*func(idz);
    end
    Xguess=Xguess./sum(Xguess(:)).*sum(blur_image(:))./(size(blur_image,3)*size(blur_image,4));
    
    a = tic;
    switch Methods
        case 'QLFM_ideal'
            Xguess = deconvGPU_QLFMwithoutScatter(maxIter,Xguess,blur_image,psf,1,UPweight);
        case 'QLFM_cali'
            Xguess = deconvGPU_QLFMwithoutScatter(maxIter,Xguess,blur_image,psf,1,UPweight);
        case 'QLFM_sca_LTR'
            Xguess = deconvGPU_QLFMScatter_L2R(maxIter,Xguess,blur_image,psf,index_1,index_2,alpha,zupdate,Vthres,UPweight);
        case 'QLFM_sca_RTL'
            Xguess = deconvGPU_QLFMScatter_R2L(maxIter,Xguess,blur_image,psf,index_1,index_2,alpha,zupdate,Vthres,UPweight);
        case 'QLFM_sca_2D'
            Xguess = deconvGPU_QLFMScatter_2D(maxIter,Xguess,blur_image,psf,index_1,index_2,alpha,zupdate,Vthres,UPweight);
        otherwise
            msg = 'No method found';
            error(msg);
    end
    
    UP_Xguess = SampleV(Xguess,index_1,index_2);
    tt = toc(a);
    disp(['frame ',num2str(frame),'multiscale reconstruction take ',num2str(tt),' sec.........']);
    
    temp=gather(UP_Xguess(CUTXY_BORDER+1:end-CUTXY_BORDER,CUTXY_BORDER+1:end-CUTXY_BORDER,CUTZ_BORDER+1:end-CUTZ_BORDER));
    imwrite3dTIFF(single(temp),[savepath,'/c1_frame',num2str(frame),'.tif']);
end



function XguessOri = deconvGPU_QLFMScatter_L2R(maxIter,Xguess,blur_image,psf,index_1,index_2,alpha,zupdate,Vthres,UPweight)
% Deconvolution with QLFM model considering about forward scattering
% Input:
    % maxIter: the maxium iteration number
    % Xguess: the estimated volume
    % blur_image: light field sub-aperture components
    % psf: light field sub-aperture point spread functions
    % index_1: downsampled index
    % index_2: upsampled index
    % alpha: learning rate
    % zupdate: update along z axis
    % Vthres: threshold for scatter coefficient
    % UPweight: update rate
% Output:
    % XguessOri: Updated 3d volume without scattering
% Rsference:

    
[img_r,img_c,~] = size(blur_image);
[psf_r,psf_c,allu,allz] = size(psf);
weight = sum(sum(sum(psf,1),2),4);
weight=weight./sum(weight(:));
weight1 = ones(13,13);
for u=1:13
    for v=1:13
        if ((u-7)^2+(v-7)^2)>=36
            weight1(u,v)=0;
        end
    end
end
ra = 0;
if(allu==45)
    ra=4;
elseif(allu==69)
    ra=5;
else
    ra=6;
end

weight2 = ones(13,13);

for u=1:13
    for v=1:13
        if ((u-7)^2+(v-7)^2)>=ra^2
            weight2(u,v)=0;
        end
    end
end



blur_image = gpuArray(single(blur_image));
XguessOri=gpuArray(single(Xguess));
XguessScatter = XguessOri;
XguessAll = XguessOri;
lambda = 525e-9;
k0 = 2*pi/lambda;
xysize = (100/13/46.52*3.25)*1e-6;
[xx,yy] = meshgrid(-(img_r-1)/2:(img_r-1)/2,-(img_r-1)/2:(img_r-1)/2);
green = gpuArray.zeros(size(XguessOri));
for idn = length(index_1):-1:2
    dz = index_1(idn)-index_1(idn-1);
    r = sqrt((xx.*xysize).^2+(yy.*xysize).^2+(dz.*0.3.*1e-6).^2);
    tmp = -exp(1j.*k0.*r.*1.4)./(4*pi.*r);
    tmp = abs(tmp).^2;
    tmp = tmp./(sum(tmp(:)));
    green(:,:,idn) = tmp;
end
VScatter = gpuArray.zeros(size(XguessAll));
Vnew = VScatter;

psf = single(psf);
largepsf = gpuArray.zeros(img_r,img_c,allz,'single');

load('map.mat');
% iteTime = tic;
flag_vup = 0;

for i=1:maxIter
    flag_uv = 0;
    for u_2=1:13
        %         bbb = tic;
        for v_2=1:13
            
            u=map_u((u_2-1)*13+v_2);
            v=map_v((u_2-1)*13+v_2);
            
            if weight2(u,v)==0
                continue;
            else
                flag_uv = flag_uv+1;
            end
            
            if weight1(u,v)==0
                continue;
            else
                %                 aa = tic;
                %                 flag_uv = flag_uv+1;
                if(i==1 && flag_uv<10)
                    tmp = gpuArray(psf(:,:,flag_uv,:));
                    largepsf((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:) ...
                        = squeeze(tmp);
                    
                    tmp = single(squeeze(blur_image(:,:,u,v)));
                    HXguess = forwardProjectACC_fft2(largepsf,XguessOri);
                    HXguessBack = backwardProjectACC_fft2(largepsf,HXguess);
                    errorBack = backwardProjectACC_fft2(largepsf,tmp);
                    
                    errorBack = (real(errorBack./HXguessBack));
                    XguessOri = XguessOri.*errorBack.*weight(flag_uv)*(1/max(weight(:))*UPweight)+(1-weight(flag_uv)*1/max(weight(:))*UPweight).*XguessOri;
                    XguessOri(isnan(XguessOri)) = 1e-5;
                    XguessOri(XguessOri<1e-5) = 1e-5;
                    XguessOri(isinf(XguessOri)) = 1e-5;
                    XguessOri = real(XguessOri);
                    XguessAll = XguessOri;
                    
                elseif(mod(flag_uv,10)==0)
                    flag_vup = flag_vup+1;
                    deltaVS = 1e100;
                    gracount  = 0;
                    while(sum(deltaVS(:))>1e-4*size(XguessAll,1)*size(XguessAll,2)*size(XguessAll,3) && gracount<=10)
                        gracount = gracount+1;
                        if(gracount==1)
                            XguessScatterNV = gpuArray.zeros(size(XguessOri));
                            %                             tmpX = XguessOri(:,:,end).*0;
                            for idn = size(XguessOri,3):-1:1
                                tmp_curr = XguessOri(:,:,idn);
                                if(idn==(size(XguessOri,3)))
                                    XguessScatterNV(:,:,idn) = tmp_curr;
                                else
                                    tmp2 = green(:,:,idn+1);
                                    tmp3 = XguessScatterNV(:,:,idn+1);
                                    tmp = forwardProjectACC_fft2(tmp2,tmp3);
                                    XguessScatterNV(:,:,idn) = tmp+tmp_curr;
                                end
                            end
                        end
                        
                        tmp = gpuArray(psf(:,:,flag_uv,:));
                        largepsf((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:) ...
                            = squeeze(tmp);
                        
                        HX = forwardProjectACC_fft2(largepsf,XguessAll);
                        Cfunc = HX-single(squeeze(blur_image(:,:,u,v)));
                        tmp= Cfunc.^2;
                        disp(['Cfunc = ',num2str(sum(tmp(:)))]);
                        dVS = gpuArray.zeros(size(VScatter));
                        for idn = size(XguessOri,3):-1:1
                            tmp1 = rot90(squeeze(largepsf(:,:,idn)),2);
                            dVS(:,:,idn) = forwardProjectACC_fft2(tmp1,Cfunc).*XguessScatterNV(:,:,idn);
                        end
                        tmp = dVS.^2;
                        disp(['gracount',num2str(gracount),'_norm gradiant = ',num2str(sum(tmp(:)))]);
                        VScatter = VScatter-alpha.*dVS;
                        VScatter(VScatter<0) = 0;
                        VScatter(VScatter>Vthres) = Vthres;
                        XguessScatter = VScatter.*XguessScatterNV;

                        XguessAll = XguessOri+XguessScatter;

                    end
                    if(flag_vup==1)
                        Vnew = VScatter;
                    else
                        VScatter = 0.1.*Vnew+VScatter;
                        Vnew = VScatter;
                    end

                    XguessAll = XguessOri+XguessScatter;
                    HX = forwardProjectACC_fft2(largepsf,XguessAll);
                    tmp = single(squeeze(blur_image(:,:,u,v)));
                    HXguessBack = backwardProjectACC_fft2(largepsf,HX);
                    errorBack = backwardProjectACC_fft2(largepsf,tmp);
                    errorBack = (real(errorBack./HXguessBack));
                    if(i==1)
                        XguessAll = XguessAll.*errorBack.*weight(flag_uv)*(1/max(weight(:))*UPweight)+...
                            (1-weight(flag_uv)*1/max(weight(:))*UPweight).*XguessAll;
                    else
                        tmp = XguessAll;
                        XguessAll = XguessAll.*errorBack.*weight(flag_uv)*(1/max(weight(:))*UPweight)+...
                            (1-weight(flag_uv)*1/max(weight(:))*UPweight).*XguessAll;
                        XguessAll = (1-zupdate.*0.5).*tmp+zupdate.*0.5.*XguessAll;
                    end
                    
                    XguessAll(isnan(XguessAll)) = 1e-5;
                    XguessAll(XguessAll<1e-5) = 1e-5;
                    XguessAll(isinf(XguessAll)) = 1e-5;
                    XguessAll = real(XguessAll);
                    XguessOri = XguessAll-XguessScatter;%有点问题的
                    XguessOri(isnan(XguessOri)) = 1e-5;
                    XguessOri(XguessOri<1e-5) = 1e-5;
                    XguessOri(isinf(XguessOri)) = 1e-5;
                    XguessOri = real(XguessOri);
                else
                    tmp = gpuArray(psf(:,:,flag_uv,:));
                    largepsf((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:) ...
                        = squeeze(tmp);
                    
                    HX = forwardProjectACC_fft2(largepsf,XguessAll);
                    tmp = single(squeeze(blur_image(:,:,u,v)));
                    HXguessBack = backwardProjectACC_fft2(largepsf,HX);
                    errorBack = backwardProjectACC_fft2(largepsf,tmp);
                    errorBack = (real(errorBack./HXguessBack));
                    if(i==1)
                        XguessAll = XguessAll.*errorBack.*weight(flag_uv)*(1/max(weight(:))*UPweight)+...
                            (1-weight(flag_uv)*1/max(weight(:))*UPweight).*XguessAll;
                    else
                        tmp = XguessAll;
                        XguessAll = XguessAll.*errorBack.*weight(flag_uv)*(1/max(weight(:))*UPweight)+...
                            (1-weight(flag_uv)*1/max(weight(:))*UPweight).*XguessAll;
                        XguessAll = (1-zupdate.*0.5).*tmp+zupdate.*0.5.*XguessAll;
                    end
                    XguessAll(isnan(XguessAll)) = 1e-5;
                    XguessAll(XguessAll<1e-5) = 1e-5;
                    XguessAll(isinf(XguessAll)) = 1e-5;
                    XguessAll = real(XguessAll);
                    XguessOri = XguessAll-XguessScatter;
                    XguessOri(isnan(XguessOri)) = 1e-5;
                    XguessOri(XguessOri<1e-5) = 1e-5;
                    XguessOri(isinf(XguessOri)) = 1e-5;
                    XguessOri = real(XguessOri);
                    XguessOri(XguessOri>1e3) = 1e3;
                    
                end
            end
            disp(['  iter ' num2str(i) ' | ' num2str(maxIter),' (u=',num2str(u), ', v=',num2str(v), '),  Energy=' num2str(sum(XguessOri(:)))]);
        end
    end
end
end

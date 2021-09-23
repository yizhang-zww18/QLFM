function Xguess = deconvGPU_QLFMwithoutScatter(maxIter,Xguess,blur_image,psf,f_shift,UPweight)
% Deconvolution with QLFM model without considering about scattering
% Input:
    % maxIter: the maxium iteration number
    % Xguess: the estimated volume
    % blur_image: light field sub-aperture components
    % psf: light field sub-aperture point spread functions
    % f_shift: compensation for fftshift
    % UPweight: update rate
% Output:
    % Xguess: Updated 3d volume
% Reference:
    

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
Xguess=gpuArray(single(Xguess));
psf = single(psf);
largepsf = gpuArray.zeros(img_r,img_c,allz,'single');
load('map.mat');

for i=1:maxIter
    flag_uv = 0;
    for u_2=1:13
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
                if(u>13)
                    continue;
                else
                    tmp = gpuArray(psf(:,:,flag_uv,:));
                    largepsf((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_c-1)/2:(img_c+1)/2+(psf_c-1)/2,:) ...
                        = squeeze(tmp);
                    FFTPSF = fftn(flip((largepsf),3));
                    HXguess =sum((fftn(Xguess).*FFTPSF),3)./allz;
                    
                    FFTPSF2 = fftn(((rot90(largepsf,2))));
                    HXguessBack = (ifftn((repmat(HXguess,[1,1,allz]).*(FFTPSF2))));
                    
                    tmp = single(squeeze(blur_image(:,:,u,v)));
                    errorBack = (ifftn((fftn(tmp,[img_r,img_c,allz]).*(FFTPSF2))));
                    
                    errorBack = (real(errorBack./HXguessBack));
                    Xguess = Xguess.*errorBack.*weight(flag_uv)*(1/max(weight(:))*UPweight)+(1-weight(flag_uv)*1/max(weight(:))*UPweight).*Xguess;

                    Xguess(isnan(Xguess)) = 1e-5;
                    Xguess(Xguess<1e-5) = 1e-5;
                    Xguess = real(Xguess);
                end
            end

            disp(['  iter ' num2str(i) ' | ' num2str(maxIter),' (u=',num2str(u), ', v=',num2str(v), '),  Energy=' num2str(sum(Xguess(:)))]);
        end
        
    end

end

if(f_shift)
    Xguess = fftshift(Xguess,1);
    Xguess = fftshift(Xguess,2);
    Xguess=gather(Xguess);
else
    Xguess=gather(Xguess);
end
end

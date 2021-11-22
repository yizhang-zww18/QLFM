clc;
clear;
cutSize = 201;
idlayer_s = 1:4:201;
idlayer2 = idlayer_s(1):idlayer_s(end);

downsample_xy = ones(1,length(idlayer2))*16;



psf = zeros(cutSize,cutSize,15,15,length(idlayer_s));
count = 0;
for idy = idlayer_s
    file = ['S:\PSF\psf_zww\20200125_genepsf_3.1746x_sim_neg400T400_dz4_15Nnum_OSR5\psf_sim_3.1746_',num2str(idy),'.mat'];
    
    count = count+1;
    load(file);
    r = size(psf_z,1);
    
    for idu = 1:15
        for idv = 1:15
            tmp = squeeze(psf_z(:,:,idu,idv));
            tmp(isnan(tmp)) = 0;
            sss = sum(tmp(:));
            if(sss>0)
                nr = round((r)/downsample_xy(idy));
                if(mod(nr,2)==0)
                    nr = nr+1;
                end
                if(downsample_xy(1)==1)
                    tmp2=tmp;
                else
                    tmp2 = imresize(tmp,[nr,nr],'nearest');
                    tmp2(tmp2<0) = 0;
                    tmp2 = tmp2./(sum(tmp2(:)))*sss;
                end
                tmp3 = zeros(r,r);
                tmp3((r+1)/2-(nr-1)/2:(r+1)/2+(nr-1)/2,(r+1)/2-(nr-1)/2:(r+1)/2+(nr-1)/2) = tmp2;
                tmp3(tmp3<0) = 0;
                tmp3(isnan(tmp3)) = 0;
                psf(:,:,idu,idv,count) = single(tmp3((r+1)/2-(cutSize-1)/2:(r+1)/2+(cutSize-1)/2,(r+1)/2-(cutSize-1)/2:(r+1)/2+(cutSize-1)/2));
            end
        end
        
    end
    disp(['idy = ',num2str(idy)]);
end

psf = single(psf);
save(['20210511_psf3x_sim_down',num2str(downsample_xy(1)),'_neg400T400_dz8','.mat'],'psf','-v7.3');
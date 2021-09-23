function projection = forwardProjectACC_fft2(psf,realspace)
% accelerated forward projection
% Input: 
    % psf: point spread function
    % realspace: 3D stack
% Output:
    % projection: projected image
% Reference:

projection=gpuArray.zeros(size(realspace,1),size(realspace,2),'single');
[rr,cc,zz] = size(realspace);
tmpa = gpuArray.zeros(2*rr-1,2*cc-1);
tmpb = gpuArray.zeros(2*rr-1,2*cc-1);
for z=1:size(realspace,3)
  %%  
    a1 = realspace(:,:,z) ;
    b1 = psf(:,:,z) ;
    tmpa((end+1)/2-(rr-1)/2:(end+1)/2+(rr-1)/2,(end+1)/2-(cc-1)/2:(end+1)/2+(cc-1)/2)=a1;
    tmpb((end+1)/2-(rr-1)/2:(end+1)/2+(rr-1)/2,(end+1)/2-(cc-1)/2:(end+1)/2+(cc-1)/2)=b1;
%     clear con1;
%     con1 = conv2(a1,b1,'same');
    con1 = ifftshift(ifft2(fft2(tmpa).*fft2(tmpb)));
%     con1 = ifft2(fft2(a1).*fft2(b1));
    projection = projection + real(con1((end+1)/2-(rr-1)/2:(end+1)/2+(rr-1)/2,(end+1)/2-(cc-1)/2:(end+1)/2+(cc-1)/2));
end
end
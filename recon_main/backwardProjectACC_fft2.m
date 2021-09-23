function Backprojection = backwardProjectACC_fft2(psf_t,projection)
% accelerated backward  projection
% Input: 
    % psf_t: point spread function
    % projection: Projected images
% Output:
    % 3D stack obtained by backward projection
% Reference:



psf_t=rot90(psf_t,2);
Backprojection=gpuArray.zeros(size(projection,1),size(projection,2),size(psf_t,3),'single');
[ra,ca]=size(projection);
b1 = gpuArray.zeros(2*ra-1,2*ca-1,'single');
a1 = gpuArray.zeros(2*ra-1,2*ca-1,'single');
for z=1:size(psf_t,3)
    a1((end+1)/2-(ra-1)/2:(end+1)/2+(ra-1)/2,(end+1)/2-(ca-1)/2:(end+1)/2+(ca-1)/2) = projection ;
    b1((end+1)/2-(ra-1)/2:(end+1)/2+(ra-1)/2,(end+1)/2-(ca-1)/2:(end+1)/2+(ca-1)/2) = psf_t(:,:,z) ;
    clear con1;
    con1 = ifftshift(ifft2(fft2(a1) .* fft2(b1)));
    Backprojection(:,:,z) = real(con1((end+1)/2-(ra-1)/2:(end+1)/2+(ra-1)/2,(end+1)/2-(ca-1)/2:(end+1)/2+(ca-1)/2));
end
end
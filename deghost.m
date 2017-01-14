function deghost()

% Dependency on Patch-GMM prior
addpath('epllcode');

% Dependency on bounded-LBFGS Optimization package
addpath('lbfgsb/lbfgsb3.0_mex1.2');

I=imread('apples.png');
I=im2double(I);
% Estimate kernel parameters

[dx dy c]=est_kernel_params(I);
padding=ceil(norm([dx dy]))+10;
[h w ch]=size(I);

I_t=I;
I_r=I;
for i=1:ch
    fprintf('Channel %d .....', i);
    %apply optimization to each channel respectively
    [it ir]=patch_gmm(I(:,:,i),h,w,c,dx, dy,i);     
    %improve color flavor by post processing
    [I_t(:,:,i), I_r(:,:,i)]=postprocess(it,ir,padding,I(:,:,i));
    
end

%Store output
imwrite(I, 'input.png');
imwrite(I_t, 't.png');
imwrite(I_r, 'r.png');
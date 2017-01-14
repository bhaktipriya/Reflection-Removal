function [ I_t I_r ] = postprocess(I_t,I_r,padding,I_in)
% Center image by subtracting with min value

I_t_main=I_t(padding+1:end-padding-1,padding+1:end-padding);
minval=min(I_t_main(:));
I_t=I_t-minval;


I_r_main=I_r(padding+1:end-padding-1,padding+1:end-padding);
minval=min(I_r_main(:));
I_r=I_r-minval;

% Match the global color flavor of the transmitted image to the original image 
sig=sqrt(sum((I_in-mean(I_in(:))).^2)/sum((I_t-mean(I_t(:))).^2));
I_t=sig*(I_t-mean(I_t(:))) + mean(I_in(:)); 

end


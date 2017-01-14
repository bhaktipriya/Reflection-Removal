function [dx,dy,c] = est_kernel_params(I)
%I is rgb image
I=rgb2gray(I);
%define a laplacian filter
filter_lp=[0 -1 0; -1 4 -1; 0 -1 0];
img_lp=imfilter(I,filter_lp);

%compute autocorrelation
auto_corr=xcorr2(img_lp);

%compute local max filter in 5x5 window
lmax1=ordfilt2(auto_corr,25, true(5));
lmax2=ordfilt2(auto_corr,24, true(5));

%Discard local maxima in neighborhoods where the first and
%second maxima are closer than a predefined threshold

% This removes local maxima caused due to
%locally flat or repetative structures.

threshold=70;
coords=(lmax1-lmax2)>=threshold; 

%We also remove local maxima within 4 pixels of origin
coords(end/2-4:end/2+4,end/2-4:end/2+4)=0;

%We select the largest one among the local maximas as 
%the Ghosting distance

lmax=lmax1.*coords;
[d_len,d]=max(lmax(:));
[dy dx]=ind2sub(size(lmax),d);
%Subtract with offset
dy=floor(size(I,1)-dy);
dx=floor(size(I,2)-dx);
dx=30;
dy=0;

%corners=corner(I);
c=0;
c=est_attenuation(I,dx,dy);


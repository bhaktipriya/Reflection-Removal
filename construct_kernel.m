function [ kernel ] = construct_kernel(h,w,dx,dy,c)

indx=[1:h*w];
indx=reshape(indx,[h,w]);

%applies shift dx dy, to show which pixels have been shifted
neigh_indx=ncircshift(indx,[dy,dx]);

%Capture indices which experience a shift and apply attenuation to them

ind=ones(h,w);
%s1 indicates pixels which have a shift and to which pixels they get
%mapped to
% the matrix dims are n*n, n=h*w
s1=sparse(indx(:),indx(:),ind);

indc=ones(h,w);
indc(neigh_indx==0)=0;
indc=indc.*c;

%s2 indicates pixels which have undergone attenuations
% the matrix dims are n*n, n=h*w
s2=sparse(indx(:),indx(:),indc);

kernel=s1+s2;

end


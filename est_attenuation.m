function [ c ] = est_attenuation( I, dx, dy )
%Find harris corners to est attenuation
corners=corner(I);
w=5;
weight=zeros(size(corners,1));
attn=zeros(size(corners,1));
for i = 1:size(corners,1)
    a=corners(i,:);
    x=a(1);
    y=a(2);
    p1=get_patch(I,x,y,w);
    nx=x+dx;
    ny=y+dy;
    p2=get_patch(I,nx,ny,w);
    if(isempty(p1)||isempty(p2))
        %disp('not empty');
        continue;
    end
    auto_corr=xcorr2(p1,p1);
%     p1v=var(p1);
%     p2v=var(p2);
     p1v=max(p1(:))-min(p1(:));
     p2v=max(p2(:))-min(p2(:));
    attn(i)=sqrt(p2v/p1v);
    if(attn(i)<1)
        score=-sum(sum(p1.*p2))/(sqrt(sum(sum(p1^2)))*sqrt(sum(sum(p2^2))));
        weight(i)=exp(-score/(2*(0.2^2)));
    end   
end

c= sum(weight.*attn)/sum(weight);
end


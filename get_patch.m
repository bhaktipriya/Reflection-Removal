function [ patch ] = get_patch( I,x,y,w)
% assuming w is odd
if ((x>w) && (y>w) &&((x+w)<size(I,2)) && ((y+w)<size(I,1)))  
    patch=I(y-w:y+w,x-w:x+w);
else
    patch=[];
end


end


function y = shift(x,axis,direction)
%SHIFT Periodic nearest-neighbor operator, shared by fluxes and projection.
s=zeros(1,ndims(x));s(axis)=direction;y=circshift(x,s);
end

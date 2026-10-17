function mask = circle_mask(coords,width)
%CIRCLE_MASK Compact smooth approximation of the exact unit-diameter cylinder.
s=sqrt(coords{1}.^2+coords{2}.^2)-.5;mask=zeros(size(s));mask(s<=-width)=1;take=abs(s)<width;mask(take)=.5-.5*sin(pi*s(take)/(2*width));
end

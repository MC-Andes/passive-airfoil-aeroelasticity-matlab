function result = channel(points,height,nu,pressureGradient,rho)
%CHANNEL Fully developed no-slip finite-volume Poiseuille verification.
% Cell-centered wall ghosts enforce u_wall=0; same central viscous stencil.
dy=height/points;y=((1:points)'-.5)*dy;L=spdiags([ones(points,1),-2*ones(points,1),ones(points,1)],-1:1,points,points)/dy^2;
L(1,1)=-3/dy^2;L(end,end)=-3/dy^2;u=L\(-pressureGradient/(rho*nu)*ones(points,1));exact=pressureGradient/(2*rho*nu)*y.*(height-y);
result=struct('y',y,'velocity',u,'exact',exact,'relativeError',norm(u-exact,inf)/max(exact), ...
    'flowRate',sum(u)*dy,'exactFlowRate',pressureGradient*height^3/(12*rho*nu));
end

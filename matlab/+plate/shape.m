function [N,dN] = shape(xi,eta)
%SHAPE Bilinear quad, counterclockwise nodes (-1,-1),(1,-1),(1,1),(-1,1).
N=.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
dN=.25*[-(1-eta),1-eta,1+eta,-(1+eta); ...
    -(1-xi),-(1+xi),1+xi,1-xi];
end

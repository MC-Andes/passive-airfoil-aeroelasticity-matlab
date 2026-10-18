function [U,p] = analytic(g,t,nu,kind)
%ANALYTIC Independent periodic Navier-Stokes solutions on native MAC locations.
d=g.dimension;U=cell(1,d);
if strcmp(kind,'taylor_green')
    assert(d==2,'nsflow:Analytic','Taylor-Green case here is two-dimensional.');
    X=g.velocityCoordinates;U{1}=sin(X{1,1}).*cos(X{1,2})*exp(-2*nu*t);
    U{2}=-cos(X{2,1}).*sin(X{2,2})*exp(-2*nu*t);
    P=g.pressureCoordinates;p=-.25*(cos(2*P{1})+cos(2*P{2}))*exp(-4*nu*t);
elseif strcmp(kind,'abc')
    assert(d==3,'nsflow:Analytic','ABC case is three-dimensional.');X=g.velocityCoordinates;
    U{1}=(sin(X{1,3})+cos(X{1,2}))*exp(-nu*t);
    U{2}=(sin(X{2,1})+cos(X{2,3}))*exp(-nu*t);
    U{3}=(sin(X{3,2})+cos(X{3,1}))*exp(-nu*t);
    P=g.pressureCoordinates;u=sin(P{3})+cos(P{2});v=sin(P{1})+cos(P{3});w=sin(P{2})+cos(P{1});
    p=-.5*(u.^2+v.^2+w.^2)*exp(-2*nu*t);p=p-mean(p,'all');
else,error('nsflow:Analytic','Unknown reference solution.');end
end

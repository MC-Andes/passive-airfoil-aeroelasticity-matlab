function R = rhs(g,U,nu,force)
%RHS Conservative centered MAC momentum fluxes and molecular viscosity.
% force{i} is acceleration on velocity-i faces. No eddy viscosity or turbulence model.
d=g.dimension;R=cell(1,d);
for i=1:d
    convection=zeros(g.count);diffusion=convection;
    for j=1:d
        diffusion=diffusion+(nsflow.shift(U{i},j,1)-2*U{i}+nsflow.shift(U{i},j,-1))/g.h(j)^2;
        if i==j
            flux=((U{i}+nsflow.shift(U{i},i,-1))/2).^2;
            convection=convection+(flux-nsflow.shift(flux,i,1))/g.h(i);
        else
            flux=.25*(U{i}+nsflow.shift(U{i},j,1)).*(U{j}+nsflow.shift(U{j},i,1));
            convection=convection+(nsflow.shift(flux,j,-1)-flux)/g.h(j);
        end
    end
    R{i}=-convection+nu*diffusion+force{i};
end
end

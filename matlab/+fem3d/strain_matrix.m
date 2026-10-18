function B = strain_matrix(grad)
%STRAIN_MATRIX Engineering shear strains, interleaved Cartesian DOFs.
B=zeros(6,3*size(grad,1));
for a=1:size(grad,1)
    g=grad(a,:);
    B(:,3*a-2:3*a)=[g(1),0,0;0,g(2),0;0,0,g(3); ...
        g(2),g(1),0;0,g(3),g(2);g(3),0,g(1)];
end
end

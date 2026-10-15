function velocity = segment_velocity(points,a,b,core)
%SEGMENT_VELOCITY Biot-Savart, unit circulation, regularized finite segment.
r1 = points-a; r2 = points-b; r0 = b-a;
n1 = vecnorm(r1,2,2); n2 = vecnorm(r2,2,2);
product = cross(r1,r2,2);
denominator = sum(product.^2,2)+core^2*sum(r0.^2);
valid = n1>core*1e-3 & n2>core*1e-3 & denominator>realmin;
velocity = zeros(size(points));
factor = sum((r1(valid,:)./n1(valid)-r2(valid,:)./n2(valid)).*r0,2);
velocity(valid,:) = product(valid,:).*(factor./denominator(valid))/(4*pi);
end

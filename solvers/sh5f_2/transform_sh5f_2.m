function [q,p,cen, Rtransform, Rtransform2, shift, Ktransform] = transform(q,p,cen)
%% Coordinate transforms Pin hole camera
a = q(1:2,1); scale = norm(a,2); b = [scale;0];
at = [-a(2);a(1)]; bt=[-b(2);b(1)];
Rtransform = [b,bt]/[a,at];
Rtransform = [Rtransform,zeros(2,1);zeros(1,2),1];
Ktransform = diag([1/scale,1/scale,1]);

q = Ktransform*Rtransform * q; q = q ./ q(3,:);
%% COORdinate transforms gen. camera
shift = cen(:,1);
a = p(:,1);    a = a/norm(a,2);
b = [0;0;1]; b = b/norm(b,2);
v = cross(a,b);    c = dot(a,b);    vx = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
Rtransform2 = eye(3) + vx + vx^2 *(1/(1+c));
cen = Rtransform2 * (cen - shift);
%% All cameras
p = Rtransform2 * p;
p = p./sqrt(sum(p.^2));
for i = 1:5
    cen(:,i) = cross(p(:,i),cen(:,i));
end
end
%% A function to transform the coordinate systems for the two cameras
function [q_d, p_d, c_d, Rtransform1, Rtransform2, shift] = transform_sh5_3(q_d, p_d, c_d)
%% Coordinate transform Q, pin hole camera
a = q_d(:,1);    a = a/norm(a,2);
b = [0;0;1]; b = b/norm(b,2);
v = cross(a,b);    c = dot(a,b);    vx = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
Rtransform1 = eye(3) + vx + vx^2 *(1/(1+c));
q_d = Rtransform1 * q_d; q_d = q_d ./ q_d(3,:);
q_d = q_d./sqrt(sum(q_d.^2));
%% Coordinate transform P, gen. camera
shift = c_d(:,1);
a = p_d(:,1);    a = a/norm(a,2);
b = [0;0;1]; b = b/norm(b,2);
v = cross(a,b);    c = dot(a,b);
vx = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
Rtransform2 = eye(3) + vx + vx^2 *(1/(1+c));
c_d = Rtransform2 * (c_d - shift);
p_d = Rtransform2 * p_d;
p_d = p_d./sqrt(sum(p_d.^2));
%% This is just a simple change in how the centers of projection are 
%% expected by the solver.
for i = 1:5
    c_d(:,i) = cross(p_d(:,i),c_d(:,i));
end
end
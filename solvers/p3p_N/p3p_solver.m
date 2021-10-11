%% Knowing N in the Gen. coordinate system, estimate the 3D points in Gen. coordinate system and execute the p3p solver.
function [Rs, ts] = p3p_solver(q, p, c, N)
%% Get the depths of the points knowing the plane and the ray and the center of projection
b = arrayfun(@(i) (-1 - N'*c(:,i))/(N'*p(:,i)), 1:3 );
% Get the 3D point, knowing the depths. 
Points_3D = [b .* p + c];

q = q ./ sqrt(sum(q.^2));
[Rs,ts] = mexp3p([q,Points_3D]);
for i=1:size(Rs,2)/3
    Rs(:, 3*i-2:3*i) = Rs(:, 3*i-2:3*i)';
    ts(:, i) = -Rs(:, 3*i-2:3*i)*ts(:, i);
end

end
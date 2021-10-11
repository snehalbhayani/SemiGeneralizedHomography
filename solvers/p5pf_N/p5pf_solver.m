%% Knowing N in the Gen. coordinate system, estimate the 3D points in Gen. coordinate system and execute the p3p solver.
function [Rs, ts, fs] = p5pf_solver(q, p, c, N)
%% Get the depths of the points knowing the plane and the ray and the center of projection
b = arrayfun(@(i) (-1 - N'*c(:,i))/(N'*p(:,i)), 1:5 );
% Get the 3D point, knowing the depths. 
Points_3D = [b .* p + c];

% R_gt= R_gt';
% t_gt = -R_gt*t_gt;

% We just need to shift the gen. coordinate system as the p4pf solver works
% by assuming that the 3rd coordinate of the input 3D points will be 0
% shift = [0;0;Points_3D(3,1)];
% Points_3D = Points_3D-shift;
% t_gt = t_gt + R_gt*shift;

q = q(1:3,:) ./ q(3,:);
% arrayfun(@(i) norm(cross((K_gt\[q(:,i);1]), R_gt*Points_3D(:,i)+t_gt),2), 1:4)

[Ps, fs] = mexp5pf([q, Points_3D]);
Rs = []; ts = [];
for i = 1:length(fs)    
    P = Ps(:,4*i-3:4*i);
    R = P(:,1:3);
    t = P(:,4);
     R = R';
     t = -R * t;
    Rs = [Rs, R];
    ts = [ts, t];
end

end
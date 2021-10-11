% Assuming 3 2d-2d corrs between two cameras in Gen. setup

function [N] = compute_N_gen(p1, p2, K_gen1, R_gen1, t_gen1, K_gen2, R_gen2, t_gen2)
p1 = (K_gen1) \ p1;
p2 = (K_gen2) \ p2;

R_12 = R_gen2*R_gen1';
t_12 = t_gen2 - R_gen2*R_gen1'*t_gen1;

c1 = cross(p2(:,1), R_12*p1(:,1)); c2 = cross(p2(:,1), t_12);
b(1) = -c2(1)/c1(1);
c1 = cross(p2(:,2), R_12*p1(:,2)); c2 = cross(p2(:,2), t_12);
b(2) = -c2(1)/c1(1);
c1 = cross(p2(:,3), R_12*p1(:,3)); c2 = cross(p2(:,3), t_12);
b(3) = -c2(1)/c1(1);

Points_3D = [b .* p1];
% Get the plane normal.
N = cross(Points_3D(:,1)- Points_3D(:,2), Points_3D(:,1)-Points_3D(:,3));
N = N/norm(N,2);
intercept = - N'*Points_3D(:,1);
N = N/intercept;

% We need to convert the plane vector into the one into the global
% coordinate system
N = (R_gen1' * N) / (1+ ( N)'*t_gen1);

end
function [Hs, Nss, Rs, ts] = sh45_3(q,p,c)
[q,p,c, Rtransform1, Rtransform2, shift]= transform_sh45_3(q,p,c);

M = nullspace_parameterization_sh45_2(q,p,c);

data = get_coeffs_sh45_3(M);
data(1,:) = data(1,:)/norm(data(1,:),2);
data(2,:) = data(2,:)/norm(data(2,:),2);
data(3,:) = data(3,:)/norm(data(3,:),2);
data(4,:) = data(4,:)/norm(data(4,:),2);
data(5,:) = data(5,:)/norm(data(5,:),2);
data = reshape(data,105,1);
sols = gb_solver_sh45_3(data);
%% Extract relative pose from the nullspace parameterization
sols = [sols;ones(1,size(sols,2))];
vecs = M * sols;
[Hs, Nss] = extract_homographies_sh45_3(vecs, sols, Rtransform1, Rtransform2, shift);
end

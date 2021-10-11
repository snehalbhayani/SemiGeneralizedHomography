function [Hs, Nss] = sh5_2(q,p,c)
[q,p,c, Rtransform1, Rtransform2, shift]= transform_sh5_2(q,p,c);
%% Computing the coefficients
M = nullspace_parameterization_sh5_2(q,p,c);
%% Companion matrix solution
data = reshape(M,16,1);
C = get_companion_matrix_sh5_2(data);
sols = eig(C);
% Removing the imaginary solutions
sols = sols(find(abs(imag(sols)) < 1e-6));
%% Extract relative pose from the nullspace parameterization
vecs = M * transpose([sols,ones(size(sols,1),1)]);
[Hs, Nss] = extract_homographies_sh5_2(vecs, sols, Rtransform1, Rtransform2, shift);
end

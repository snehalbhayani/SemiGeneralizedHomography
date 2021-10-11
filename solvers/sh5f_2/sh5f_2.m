function [Hs, fss, Nss] = sh5f_2(q,p,c)
[q,p,c, Rtransform, Rtransform2, shift, Ktransform] = transform_sh5f_2(q,p,c);

%% COmputing the coefficients
M = nullspace_parameterization_sh5f_2(q,p,c);

%% Companion matrix solution
data = reshape(M,16,1);
C = get_companion_matrix_sh5f_2(data);
sols = eig(C);
% Removing the imaginary solutions
sols = sols(find(abs(imag(sols)) < 1e-6));

%% Extract relative pose from the nullspace parameterization
vecs = M * transpose([sols,ones(size(sols,1),1)]);
[Hs, fss, Nss] = extract_homographies_sh5f_2(vecs, sols, Ktransform, Rtransform, Rtransform2, shift);
end

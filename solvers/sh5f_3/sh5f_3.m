function [Hs, fs, Ns] = sh5f_3(q,p,c)
%% Transform the coordinate system
[q,p,c, Rtransform1, Rtransform2, shift, Ktransform] = transform_sh5f_3(q,p,c);
%% Obtain the nullspace parameterization
M = nullspace_parameterization_sh5f_3(q,p,c);
h = -M(:,1:8)\M(:,10);
h = h(1:6);
b = M([6,8], 7:9);
A = -M([6,8], [1:6,10]) * [h;1];
c1 = b(:,1:2)\A; c2 = -b(:,1:2)\b(:,3);
data = [c1;c2;h];
C = get_companion_matrix_sh5f_3(data);
sols = eig(C);

%% Removing the imaginary solutions
% in the worst case we will have 3 real solutions
sols = sols(find(abs(imag(sols)) < 1e-6));
vecs = zeros(8,length(sols));
for i=1:length(sols)
    s = sols(i);
    v = data(1:2) + data(3:4) * s;
    vecs(:,i) = [h;v];
end
%% Extract the semi-generalized homographies
[Hs, fs, Ns] = extract_homographies_sh5f_3(vecs, sols, Ktransform, Rtransform1, Rtransform2, shift);

end

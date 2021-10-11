function [Hs, fss, Nss] = extract_homographies_sh5f_3(vecs, sols, Ktransform, Rtransform, Rtransform2, shift)
vecs = [vecs;sols'];
Hs = [];
fss = [];
Nss = [];
for i = 1:size(vecs,2)
    vec = vecs(:,i);
    [h11, h21, h12, h22, h32, h33, nx, ny, nz] = feval(@(x) x{:}, num2cell(vec));
    Nscaled = [nx;ny;nz];
    Gscaled = [h11,h12,-h11;h21,h22,-h21;1,h32,h33];
    scale = sqrt((nx ^ 2 + ny ^ 2) / (h11 ^ 2 * ny ^ 2 - 2 * h11 * h12 * nx * ny + h12 ^ 2 * nx ^ 2 + h21 ^ 2 * ny ^ 2 - 2 * h21 * h22 * nx * ny + h22 ^ 2 * nx ^ 2 + h32 ^ 2 * nx ^ 2 - 2 * h32 * nx * ny + ny ^ 2));
    % We need to ensure that the plane vector has a sign such that the
    % depth of the point is +ve.
    % If the 3D point is X = alpha * K^-1 * q then,
    % alpha * N^T * K^-1 * q + 1 = 0
    % which means that N^T * K^-1 * q must be -ve so that alpha is +ve
    scale = -sign(scale*Nscaled'*[1;0;1]) * scale;
    Ns = scale * Nscaled;
    Gs = scale * Gscaled;
    
    % Extracting the focal length.
    ws = abs(ny/sqrt(scale^2*h11^2*ny^2 + 2*h11*h12*ny*nz*scale^2 + h12^2*nz^2*scale^2 + h21^2*ny^2*scale^2 + 2*h21*h22*ny*nz*scale^2 ...
        + h22^2*nz^2*scale^2 + h32^2*nz^2*scale^2 -2*h32*h33*ny*nz*scale^2 + h33^2*ny^2*scale^2 - nz^2));
    K = diag([1,1,ws]);
    
    H = Gs * K;
    N = K * Ns;
    ws = 1/ws;
    
    % Homography solution is reverse transformed
    ws = ws/Ktransform(1,1);
    Ns = Rtransform' * N;
    H = Rtransform2' * H *  Rtransform  - shift * transpose(Ns);
    fss = [fss, ws];
    Hs = [Hs, H];
    Nss = [Nss, Ns];
    
end

end

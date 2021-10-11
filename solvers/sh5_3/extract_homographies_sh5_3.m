function [Hs, Nss] = extract_homographies_sh5_3(vecs, sols, Rtransform, Rtransform2, shift)
vecs = [vecs;sols'];
Nss = [];
Hs = [];

for i = 1:size(vecs,2)
    vec = vecs(:,i);
    [h11, h21, h12, h22, h32, h33, nx, ny, nz] = feval(@(x) x{:}, num2cell(vec));
    
    Nscaled = [nx;ny;nz];
    Gscaled = [h11,h12,0;h21,h22,0;1,h32,h33];
    % Using the elimination ideal(the simplest polynomial before we eliminated h31 from it) to extract h31 from it.
    scale = sqrt((ny ^ 2 + nz ^ 2) / (h12 ^ 2 * nz ^ 2 + h22 ^ 2 * nz ^ 2 + h32 ^ 2 * nz ^ 2 - 2 * ny * h32 * h33 * nz + ny ^ 2 * h33 ^ 2));
    
    % We need to ensure that the plane vector has a sign such that the
    % depth of the point is +ve.
    % If the 3D point is X = alpha * q then,
    % alpha * N^T * q + 1 = 0
    % which means that N^T * q must be -ve so that alpha is +ve
    scale = -sign(scale*Nscaled'*[0;0;1]) * scale;
    Ns = scale * Nscaled;
    H = scale * Gscaled;
    
    % The reverse coordinate transform for
    % homography and the plane vector
    Ns = Rtransform' * Ns;    
    H = Rtransform2' * H *  Rtransform  - shift * transpose(Ns);

    Hs = [Hs, H];      
    Nss = [Nss, Ns];
        
end
end


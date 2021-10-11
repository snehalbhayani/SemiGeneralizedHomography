function [Hs, Nss] = extract_homographies_sh45_2(vecs, sols, Rtransform, Rtransform2, shift)
vecs = [vecs;sols];
Nss = [];
Hs = [];

% residues = zeros(size(vecs,2),1);
for i = 1:size(vecs,2)
    vec = vecs(:,i);
    [g11, g21, g12, g22, g31, g32, nx, ny, nz] = feval(@(x) x{:}, num2cell(vec));
    
    %     residues(i) = h11 ^ 2 * ny ^ 2 * nz + h11 ^ 2 * nz ^ 3 - h12 ^ 2 * nx ^ 2 * nz - h12 ^ 2 * nz ^ 3 + h21 ^ 2 * ny ^ 2 * nz + h21 ^ 2 * nz ^ 3 - h22 ^ 2 * nx ^ 2 * nz - h22 ^ 2 * nz ^ 3 - h32 ^ 2 * nx ^ 2 * nz - h32 ^ 2 * nz ^ 3 + 2 * h32 * h33 * nx ^ 2 * ny + 2 * h32 * h33 * ny * nz ^ 2 + h33 ^ 2 * nx ^ 2 * nz - h33 ^ 2 * ny ^ 2 * nz - 2 * h33 * nx * ny ^ 2 - 2 * h33 * nx * nz ^ 2 + ny ^ 2 * nz + nz ^ 3;
    % end
    % vec = vecs(:,find(abs(residues) == min(abs(residues)), 1));
    % [h11, h21, h12, h22, h32, h33, nx, ny, nz] = feval(@(x) x{:}, num2cell(vec));
    
    Nscaled = [nx;ny;nz];
    Gscaled = [g11,g12,0;g21,g22,0;g31,g32,1];
    % Using the elimination ideal(the simplest polynomial before we eliminated h31 from it) to extract h31 from it.
    scale = sqrt((ny ^ 2 + nz ^ 2) / (g12 ^ 2 * nz ^ 2 + g22 ^ 2 * nz ^ 2 + g32 ^ 2 * nz ^ 2 - 2 * ny * nz * g32 + ny ^ 2));
%     scale = sqrt((ny ^ 2 + nz ^ 2) / (h12 ^ 2 * nz ^ 2 + h22 ^ 2 * nz ^ 2 + h32 ^ 2 * nz ^ 2 - 2 * ny * h32 * h33 * nz + ny ^ 2 * h33 ^ 2));

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


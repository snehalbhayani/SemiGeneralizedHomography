function HCam = compute_H_cam_from_H_gen(HGen, K, R, t, N)
% Center of projection
center = -R'*t;
HCam = K * R * (HGen + center * N');
end
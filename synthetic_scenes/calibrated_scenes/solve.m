clear all;
clc;
addpath ../thirdparty
addpath ../tools
addpath ../common
addpath ../../solvers/sh5_3/;
addpath ../../solvers/sh5_2/;
addpath ../../solvers/sh4.5_2/;
addpath ../../solvers/sh4.5_3/;
addpath ../../solvers/p3p_N/;
addpath ../../solvers/sh5_4/;
addpath Noise_boxplot;

%% Test scenes on all solvers for different pixel noise levels
errors = {};
noise_levels = [0, 0.1,1.0,2.0];
for j = 1:size(noise_levels,2)
    for i = 1:100
        i
        all_solver_errors = test(noise_levels(j));
        tn_err(i,:) = all_solver_errors(1,:);
        rot_err(i,:) = all_solver_errors(2,:);
    end
    for k = 1:6
        errors{k}(j).tn_err = log10(tn_err(:,k));
        errors{k}(j).rot_err = (rot_err(:,k));
    end
end
%% Box Plots
boxplot_noise_rot;
boxplot_noise_tn;

%% Test each scene
function errors = test(noise_amp)

f = 5;
K_pinhole = f*diag([1,1,1/f]);
f_gen = 5;
K_gen = diag([f_gen, f_gen, 1]);

%noise
pixel = 1/1000;
noise = noise_amp * pixel;

% How many points to sample
Npoints = 1000;
Ncams = 7;

%% Various configurations of scene + motion
indicesOfPts = [1,2,3,4,5,6];

sceneType = {'randomplane' 'random'};
[Pgt M m m2] = GenerateScene(Npoints, 10, Ncams, 20, 35 , 0, noise,...
[K_pinhole;K_gen;K_gen;K_gen;K_gen;K_gen;K_gen], sceneType, [], [], [], true, 0);
for i = 1:Ncams
    temp = Pgt{i} * [M; ones(1,size(M,2))];
    indices = find(temp(3,:) > 0);
    M = M(:, indices);
end

%% Pinhole camera
% For all solvers we need 6 2d-2d point matches.
center1 = null(Pgt{1}); center1 = center1(1:3)/center1(4,1);
P = Pgt{1};
t = center1;
R = (K_pinhole\P(:,1:3))';

N = R'*M-R'*t;
q_d = K_pinhole\[m{1}(:,[indicesOfPts(1:6)]);ones(1,6)];


%% Generalized camera
c1 = null(Pgt{2}); c1 = c1(1:3)/c1(4);
c2 = null(Pgt{3}); c2 = c2(1:3)/c2(4);
c3 = null(Pgt{4}); c3 = c3(1:3)/c3(4);
c4 = null(Pgt{5}); c4 = c4(1:3)/c4(4);
c5 = null(Pgt{6}); c5 = c5(1:3)/c5(4);
c6 = null(Pgt{7}); c6 = c6(1:3)/c6(4);

R_gen1 = K_gen\Pgt{2}(:,1:3);
R_gen2 = K_gen\Pgt{3}(:,1:3);
R_gen3 = K_gen\Pgt{4}(:,1:3);
R_gen4 = K_gen\Pgt{5}(:,1:3);
R_gen5 = K_gen\Pgt{6}(:,1:3);
R_gen6 = K_gen\Pgt{7}(:,1:3);

R_gt = R; t_gt = t;
N_gt = cross(N(:,indicesOfPts(1))-N(:,indicesOfPts(2)), N(:,indicesOfPts(1))-N(:,indicesOfPts(3)));
N_gt = N_gt/(-N_gt'*N(:,indicesOfPts(1)));

%% Test sH5_2
p_d = [(R_gen1'/K_gen)*[m{2}(:,indicesOfPts(1));1], ...
    (R_gen2'/K_gen)*[m{3}(:,indicesOfPts(2));1],...
    (R_gen3'/K_gen)*[m{4}(:,indicesOfPts(3));1],...
    (R_gen4'/K_gen)*[m{5}(:,indicesOfPts(4));1],...
    (R_gen5'/K_gen)*[m{6}(:,indicesOfPts(5));1]];
p_d = p_d./p_d(3,:);

Hss=[]; Nss=[];
[Hss, Nss] = sh5_2([q_d(:,1),q_d(:,2),q_d(:,3),q_d(:,4),q_d(:,5)],...
    [p_d(:,1),p_d(:,2),p_d(:,3),p_d(:,4),p_d(:,5)],...
    [c1,c2,c3,c4,c5]);
Rs = []; ts = [];
for i=1:size(Nss,2)
    sols = decomp_homo(Hss(:, 3*i-2:3*i));
    t = [];
    R = [];
    N = [];
    for k = 1:length(sols)
        temp = sols(k).T;
        t = [t, temp(1:3,1:3)'*temp(1:3,4)/norm(Nss(:,i),2)];
        R = [R, temp(1:3,1:3)'];
        N = [N, sols(k).n*norm(Nss(:,i),2)];
    end
    
    for j = 1:size(t,2)
        Rs = [Rs, R(:,3*j-2:3*j)];
        ts = [ts, t(:,j)];
    end
end
ind = find( abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2))) == min(abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2)))), 1);
error_sh5_2  = [norm(ts(:,ind)-t_gt,2)/norm(t_gt,2), abs(acosd((trace(R_gt'*Rs(:,3*ind-2:3*ind))-1)/2))]';

%% Test sH5_3

p_d = [(R_gen1'/K_gen)*[m{2}(:,indicesOfPts(1));1], ...
    (R_gen1'/K_gen)*[m{2}(:,indicesOfPts(2));1],...
    (R_gen1'/K_gen)*[m{2}(:,indicesOfPts(3));1],...
    (R_gen2'/K_gen)*[m{3}(:,indicesOfPts(4));1],...
    (R_gen2'/K_gen)*[m{3}(:,indicesOfPts(5));1]];
p_d = p_d./p_d(3,:);

[Hss, Nss] = sh5_3([q_d(:,1),q_d(:,2),q_d(:,3),q_d(:,4),q_d(:,5)],...
    [p_d(:,1),p_d(:,2),p_d(:,3),p_d(:,4),p_d(:,5)],...
    [c1,c1,c1,c2,c2]);

Rs = []; ts = [];
Ns = [];
for i=1:size(Nss,2)
    sols = decomp_homo(Hss(:, 3*i-2:3*i));
    t = [];
    R = [];
    N = [];
    for k = 1:length(sols)
        temp = sols(k).T;
        t = [t, temp(1:3,1:3)'*temp(1:3,4)/norm(Nss(:,i),2)];
        R = [R, temp(1:3,1:3)'];
        N = [N, sols(k).n*norm(Nss(:,i),2)];
    end
    
    for j = 1:size(t,2)
        Rs = [Rs, R(:,3*j-2:3*j)];
        ts = [ts, t(:,j)];
        Ns = [Ns, Nss(:,i)];
    end
end
ind = find( abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2))) == min(abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2)))), 1);
error_sh5_3  = [norm(ts(:,ind)-t_gt,2)/norm(t_gt,2), abs(acosd((trace(R_gt'*Rs(:,3*ind-2:3*ind))-1)/2)) ]';

%% Test sH5_4
p_d = [(R_gen1'/K_gen)*[m{2}(:,indicesOfPts(1));1], ...
    (R_gen1'/K_gen)*[m{2}(:,indicesOfPts(2));1],...
    (R_gen1'/K_gen)*[m{2}(:,indicesOfPts(3));1],...
    (R_gen1'/K_gen)*[m{2}(:,indicesOfPts(4));1],...
    (R_gen2'/K_gen)*[m{3}(:,indicesOfPts(5));1]];
p_d = p_d./p_d(3,:);

Hss = sh5_4(q_d(:,1:5), p_d, c1, c2);

Rs = []; ts = [];
Ns = [];
i=1;
sols = decomp_homo(Hss(:, 3*i-2:3*i));
t = [];
R = [];
N = [];
for k = 1:length(sols)
    temp = sols(k).T;
    scale = scale_from_Rt(q_d(:,5), p_d(:,5),c2-c1, temp(1:3,1:3)', temp(1:3,1:3)'*temp(1:3,4), sols(k).n);
    t = [t, temp(1:3,1:3)'*temp(1:3,4)/scale];
    R = [R, temp(1:3,1:3)'];
    N = [N, sols(k).n*scale];
end

for j = 1:size(t,2)
    Rs = [Rs, R(:,3*j-2:3*j)];
    ts = [ts, t(:,j)];
    Ns = [Ns, N(:,j)];
end

ts = ts + c1;
ind = find( abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2))) == min(abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2)))), 1);
error_sh5_4  = [norm(ts(:,ind)-t_gt,2)/norm(t_gt,2), abs(acosd((trace(R_gt'*Rs(:,3*ind-2:3*ind))-1)/2)) ]';

%% Test sh4.5_2
p_d = [(R_gen1'/K_gen)*[m{2}(:,indicesOfPts(1));1], ...
    (R_gen2'/K_gen)*[m{3}(:,indicesOfPts(2));1],...
    (R_gen3'/K_gen)*[m{4}(:,indicesOfPts(3));1],...
    (R_gen4'/K_gen)*[m{5}(:,indicesOfPts(4));1],...
    (R_gen5'/K_gen)*[m{6}(:,indicesOfPts(5));1]];
p_d = p_d./p_d(3,:);

[Hss, Nss] = sh45_2([q_d(:,1),q_d(:,2),q_d(:,3),q_d(:,4),q_d(:,5)],...
    [p_d(:,1),p_d(:,2),p_d(:,3),p_d(:,4),p_d(:,5)],...
    [c1,c2,c3,c4,c5]);


Rs = []; ts = [];
for i=1:size(Nss,2)
    try
        sols = decomp_homo(Hss(:, 3*i-2:3*i));
        t = [];
        R = [];
        N = [];
        for k = 1:length(sols)
            temp = sols(k).T;
            t = [t, temp(1:3,1:3)'*temp(1:3,4)/norm(Nss(:,i),2)];
            R = [R, temp(1:3,1:3)'];
            N = [N, sols(k).n*norm(Nss(:,i),2)];
        end
        for j = 1:size(t,2)
            Rs = [Rs, R(:,3*j-2:3*j)];
            ts = [ts, t(:,j)];
            Ns = [Ns, N(:,j)];
            
        end
    catch
        1;
    end
end
try
    ind = find( abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2))) == min(abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2)))), 1);    
    error_sh45_2  = [norm(ts(:,ind)-t_gt,2)/norm(t_gt,2), abs(acosd((trace(R_gt'*Rs(:,3*ind-2:3*ind))-1)/2))]';
catch
    error_sh45_2 = [inf;inf];
end

%% Test sH4.5_3vars
p_d = [(R_gen1'/K_gen)*[m{2}(:,indicesOfPts(1));1], ...
    (R_gen1'/K_gen)*[m{2}(:,indicesOfPts(2));1],...
    (R_gen1'/K_gen)*[m{2}(:,indicesOfPts(3));1],...
    (R_gen2'/K_gen)*[m{3}(:,indicesOfPts(4));1],...
    (R_gen3'/K_gen)*[m{4}(:,indicesOfPts(5));1]];
p_d = p_d./p_d(3,:);
%
[Hss, Nss] = sh45_3([q_d(:,1),q_d(:,2),q_d(:,3),q_d(:,4),q_d(:,5)],...
    [p_d(:,1),p_d(:,2),p_d(:,3),p_d(:,4),p_d(:,5)],...
    [c1,c1,c1,c2,c3]);

    
Rs = []; ts = [];
for i=1:size(Nss,2)
    try
    sols = decomp_homo(Hss(:, 3*i-2:3*i));
    t = [];
    R = [];
    N = [];
    for k = 1:length(sols)
        temp = sols(k).T;
        t = [t, temp(1:3,1:3)'*temp(1:3,4)/norm(Nss(:,i),2)];
        R = [R, temp(1:3,1:3)'];
        N = [N, sols(k).n*norm(Nss(:,i),2)];
    end
    for j = 1:size(t,2)
        Rs = [Rs, R(:,3*j-2:3*j)];
        ts = [ts, t(:,j)];
        Ns = [Ns, N(:,j)];           
    end
    catch
    end
end
    
try
    ind = find( abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2))) == min(abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2)))), 1); 
    error_sh45_3  = [norm(ts(:,ind)-t_gt,2)/norm(t_gt,2), abs(acosd((trace(R_gt'*Rs(:,3*ind-2:3*ind))-1)/2))]';
catch
    error_sh45_3 = [inf;inf];
end
%% Test p3pN
t_gen1 = K_gen\Pgt{2}(:,4);
t_gen2 = K_gen\Pgt{3}(:,4);

p1 = [m{2}(:,indicesOfPts(1:3));ones(1,3)];
p2 = [m{3}(:,indicesOfPts(1:3));ones(1,3)];
p_d = [(R_gen4')*[m{5}(:,indicesOfPts(4));1], ...
    (R_gen5')*[m{6}(:,indicesOfPts(5));1], ...
    (R_gen6')*[m{7}(:,indicesOfPts(6));1]];

c4 = null(Pgt{5}); c4 = c4(1:3)/c4(4);
c5 = null(Pgt{6}); c5 = c5(1:3)/c5(4);
c6 = null(Pgt{7}); c6 = c6(1:3)/c6(4);

p_d = p_d./p_d(3,:);
c_d = [c4,c5,c6];

N = compute_N_gen(p1, p2, K_gen, R_gen1, t_gen1, K_gen, R_gen2, t_gen2);
[Rs, ts] = p3p_solver(q_d(:,4:6), p_d, c_d, N);


ind = find( abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2))) == min(abs(arrayfun(@(i) norm([Rs(:,3*i-2:3*i), ts(:,i)]-[R_gt,t_gt],2), 1:size(ts,2)))), 1);
try
    error_p3pn  = [norm(ts(:,ind)-t_gt,2)/norm(t_gt,2), abs(acosd((trace(R_gt'*Rs(:,3*ind-2:3*ind))-1)/2))]';
catch
    error_p3pn  = [inf;inf];
end
%% Compiling all the solvers' errors
errors = [error_sh5_2,error_sh45_2, error_sh5_3, error_sh45_3, error_sh5_4, error_p3pn];
end


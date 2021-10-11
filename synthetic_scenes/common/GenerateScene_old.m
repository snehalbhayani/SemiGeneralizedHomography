% 16.nov.2008 - pridany noise pre kazdu kameru separe
% 20.5.2009 - rozne radial distortion modely
% last update 22.oct.2008
%
%
% synthetic scene generator
% by Martin Bujnak, Jan2006

% function [P M m mnoiseless npoints] = GenerateScene(npoints, radius, ncams, minz, maxz, noutliers, anoise, camcalib, sceneType, pointsseed, camseed, rdcoefs, removeBehind)
%
% generate 3D strcture, cameras and 2D reprojections with noise and ouliers
%
% npoints - number of 3D points
% radius - size of the scene
% ncams - number of cameras
% minz - minimal distance from the scene
% maxz - maximal distance from the scene
% noutliers - number of ouliers in 0..1, 0.3 means 30%
% anoise - noise amplitude - pixels in image coordinate system
% camcalib - camera calibration matrices. Two possibilities:
%             one calib. matrix - all cameras will have the same calib
%             #ncams calib.matrices (3*N x 3)- P is modified with correspnding
%             calib. matric
% pointsseed (optional) - random generator seed for 3D points generator
% camseed (optional) - random generator seed for camera generator
% sceneType - either cell array or sigle string (motionType)
%   sceneType{1} = modelType - 3D data generator type : randombox
%   (default), regularbox
%   motionType{1} = camera configuration - random, sideway, forward, cake,
%                                          parallel axes. 
%                   Optionally a parameter may appear in () behind camear motion type
%                   string e.g. cake(10) means camera form 10' of a circle
%
% rdcoefs - radial distortion - apply radial distortion (set of coefficients for each camera, of one if all are the same)
% removeBehind - set 1 to remove points that are behind at least one camera
%                scene during generation

function [P M m mnoiseless npoints] = GenerateScene(npoints, radius, ncams, minz, maxz, noutliers, anoise, camcalib, sceneType, pointsseed, camseed, rdcoefs, removeBehind)

    screenSize = [-0.5 0.5
                  -0.5 0.5];

    if iscell(sceneType)

        modelType = sceneType{1};
        motionType = sceneType{2};
        
    else

        modelType = sceneType;
        motionType = 'random';
    end

    bAffine = false;
    if (ncams < 0)
        
        bAffine = true;
        ncams = -ncams;
    end

    if (nargin > 9) && ~isempty(pointsseed)
        
        randn('state',pointsseed);
    end

    if (size(radius, 2) < 2) 
        radius(2) = radius(1);
    end

    if (size(radius, 2) < 3)
        radius(3) = radius(1);
    end
    
    if nargin < 9 

        modelType = 'randombox';
    end
    
    modelSize = 0;
    if strncmp(modelType, 'regularbox', 10)
    
        modelSize = 10;
        side = ceil(npoints^(1/3));
        npoints = side^3;
        
        c = 1;
        for i=1:side
            for j=1:side
                for k=1:side
                    
                    M(:,c) = [(radius(1) * ((i-1)/(side-1)) - radius(1)/2);
                              (radius(1) * ((j-1)/(side-1)) - radius(1)/2);
                              (radius(1) * ((k-1)/(side-1)) - radius(1)/2)];
                          
                    c = c + 1;
                end
            end
        end
        
    elseif strncmp(modelType, 'randombox', 9)

        %generate 3D points
        modelSize = 9;
        M = radius(1)*(rand(3, npoints) - 0.5);

    elseif strncmp(modelType, 'randomplane', 11)

        %generate 3D points
        modelSize = 11;
        M = radius(1)*(rand(2, npoints) - 0.5);
        M(3, :) = 0;

    elseif strncmp(modelType, 'randomplaneX', 12)

        %generate 3D points
        modelSize = 12;
        M = radius(1)*(randn(3, npoints) - 0.5);
        M(1, :) = 0;

    elseif strncmp(modelType, 'randomplaneY',12)

        %generate 3D points
        modelSize = 12;
        M = radius(1)*(randn(3, npoints) - 0.5);
        M(2, :) = 0;

    elseif strncmp(modelType, 'regularplane', 12)

        modelSize = 12;
        side = ceil(npoints^(1/2));
        npoints = side^2;
        
        c = 1;
        for i=1:side
            for j=1:side

               M(:,c) = [(2 * radius(1) * (i/side) - radius(1));
                         (2 * radius(1) * (j/side) - radius(1))];
               c = c + 1;
            end
        end
        M(3, :) = 0;
        
    elseif strncmp(modelType, 'regularplaneX', 13)

        modelSize = 13;
        side = ceil(npoints^(1/2));
        npoints = side^2;
        
        c = 1;
        for i=1:side
            for j=1:side

               M([2 3],c) = [(2 * radius(1) * (i/side) - radius(1));
                         (2 * radius(1) * (j/side) - radius(1))];
               c = c + 1;
            end
        end
        M(1, :) = 0;
        
    elseif strncmp(modelType, 'regularplaneY', 13)

        modelSize = 13;
        side = ceil(npoints^(1/2));
        npoints = side^2;
        
        c = 1;
        for i=1:side
            for j=1:side

               M([1 3],c) = [(2 * radius(1) * (i/side) - radius(1));
                         (2 * radius(1) * (j/side) - radius(1))];
               c = c + 1;
            end
        end
        M(2, :) = 0;
    end    
    
    if length(modelType) > modelSize

        % (l%) expected
        parr = strfind(modelType, ')');
        num = str2double(modelType((modelSize+2):(parr-2)));
        
        % l% of the scene will be general 3D structure

        pts3d = npoints*(num/100);
        M3d = radius(1)*(randn(3, pts3d) - 0.5);
        M(:, (npoints-pts3d+1):npoints) = M3d;
    end
    
   
    if nargin > 10 && ~isempty(camseed)
        randn('state',camseed);
    end
    
    %create N - random cameras...
    P={};
    
    % test (l%)
    parl = strfind(motionType, '(');
    if ~isempty(parl)

        parr = strfind(motionType, ')');
        motionParam = str2double(motionType((parl+1):(parr-1)));
    else 
        motionParam = [];
    end
    
    if strncmp(motionType, 'random', 6)
    
        for i = 1:ncams

            dir = rand(3,1);
            dir = dir / norm(dir);
            pos = dir * ((rand * (maxz - minz)) + minz);
            at = radius(2) * (rand(3, 1) - 0.5);

            upDir = [rand rand rand];
            upDir = upDir / norm(upDir);
            
            
            P{i} = SetLookAt(pos', at', upDir);
            P{i}(:,1:3) = P{i}(:,1:3)';

            if bAffine

                P{i}(3,:) = [0 0 0 1];
            end
        end
        
    elseif strncmp(motionType, 'xyplane', 6)

        dir = [0;0;1];
        for i = 1:ncams

            pos = dir * ((randn * (maxz - minz)) + minz);
            at = radius(2) * (randn(3, 1) - 0.5);
            P{i} = SetLookAt(pos', at', [0 1 0]);
            P{i}(:,1:3) = P{i}(:,1:3)';

            if bAffine

                P{i}(3,:) = [0 0 0 1];
            end
        end        
        
     elseif strncmp(motionType, 'sideway', 7)

        dir = [0;0;1];
        movdir = [0;1;0];
        pos = dir * ((randn * (maxz - minz)) + minz);
        at = radius(2) * (randn(3, 1) - 0.5);
        for i = 1:ncams

            dPos = movdir * (randn * radius(3));
            pos = pos + dPos;
            at = at + dPos;
            P{i} = SetLookAt(pos', at', [0 1 0]);
            P{i}(:,1:3) = P{i}(:,1:3)';

            if bAffine

                P{i}(3,:) = [0 0 0 1];
            end
        end        
       
     elseif strncmp(motionType, 'forward', 7)

        dir = [0;0;1];
        for i = 1:ncams

            pos = dir * ((randn * (maxz - minz)) + minz);
            at = radius(2) * (randn(3, 1) - 0.5);
            P{i} = SetLookAt(pos', at', [0 1 0]);
            P{i}(:,1:3) = P{i}(:,1:3)';

            if bAffine

                P{i}(3,:) = [0 0 0 1];
            end
        end        
       
     elseif strncmp(motionType, 'parallel axes', 7)

        dir = [0;0;1];
        for i = 1:ncams

            pos = dir * ((randn * (maxz - minz)) + minz);
            at = radius(2) * (randn(3, 1) - 0.5);
            P{i} = SetLookAt(pos', at', [0 1 0]);
            P{i}(:,1:3) = P{i}(:,1:3)';

            if bAffine

                P{i}(3,:) = [0 0 0 1];
            end
        end        
       
    elseif strncmp(motionType, 'cake', 4)

        if isempty(motionParam)
            motionParam = 360;
        end

        cradius = ((randn * (maxz - minz)) + minz);
        zpos = ((randn * (maxz - minz)) + minz);
        pos = [0;0;zpos];
        for i = 1:ncams

            pos(1) = cradius*cos(pi*((i-1) * motionParam/(ncams))/180);
            pos(2) = cradius*sin(pi*((i-1) * motionParam/(ncams))/180);
            at = radius(2) * (randn(3, 1) - 0.5);
            
            P{i} = SetLookAt(pos', at', [0 1 0]);
            P{i}(:,1:3) = P{i}(:,1:3)';

            if bAffine

                P{i}(3,:) = [0 0 0 1];
            end
        end        
    end
    
    %un-calibrate cameras
    if size(camcalib, 1) == 0
        camcalib = repmat(eye(3), ncams, 1);
    elseif size(camcalib, 1) == 3
        camcalib = repmat(camcalib(1:3, 1:3), ncams, 1);
    end

    if size(camcalib, 2) > 1
        for i=1:ncams
            P{i} = camcalib(((i-1)*3+1):(i*3) ,1:3)*P{i};
        end
    end
    
    %create 2D observations...
    front=zeros(1, npoints);
    for j=1:ncams
        
        jpoints=[];
        for i = 1:npoints
        
            p = P{j} * [M(:,i);1];
            
            if p(3) < 0
                front(i) = -1;
            end
            
            p2d = p / p(3);
            jp3d(:,i) = p;
            jpoints(:,i) = p2d(1:2);
        end
        m{j} = jpoints;
        m3d{j} = jp3d;
    end
    
    if nargin > 13 && removeBehind
        
        frontal = find(front == 0);
        for j=1:ncams
            m{j} = m{j}(:, frontal);
            m3d{j} = m3d{j}(:, frontal);
        end
    end

    % save gt
    if ncams > 0
        mnoiseless = m;
    end
    
    % apply radial distrion
    if nargin > 11 && ~isempty(rdcoefs)
        
        if iscell(rdcoefs)
            
            modeltype = rdcoefs{1};
            rdcoefs = rdcoefs{2};
        else
            
            modeltype = 0;
        end

        if (modeltype == 0)
            
            % [ x ]   [     0       ]
            % [ y ] + [     0       ]
            % [ 1 ]   [ l*(x^2+y^2) ]

            % R = 1 ./ (1+lambda(i)*sum(p.^2));
            % p = p .* repmat(R, 2, 1);
            %
            %
            % eq = [x;y] .* repmat(R, 2, 1) - [xu;yu] = 0
            % res=solve([char(eq(1)) '=0'], [char(eq(2)) '=0'], 'x','y');

            % apply inverse to the div model
            for cam=1:ncams

                if length(rdcoefs) == ncams

                    lambda = rdcoefs(cam);
                else
                    lambda = rdcoefs(1);
                end       

                if lambda == 0 
                    continue;
                end

                for j=1:npoints

                    xu = m{cam}(1,j);
                    yu = m{cam}(2,j);
                    m{cam}(1,j) = 1/2*xu/(lambda*yu^2+xu^2*lambda)*(1-(1-4*lambda*yu^2-4*xu^2*lambda)^(1/2));
                    m{cam}(2,j) = 1/2/(lambda*yu^2+xu^2*lambda)*(1-(1-4*lambda*yu^2-4*xu^2*lambda)^(1/2))*yu;
                end
            end
            
        elseif (modeltype == 5) 
            
            % polynomial model - higher coefs first!
            % (x,y) = (x,y) + r(x,y)
            
            for cam=1:ncams

                if size(rdcoefs, 1) == ncams

                    coefs = rdcoefs(cam, :);
                else
                    coefs = rdcoefs(:);
                end     
                
                coefs = [coefs 1];

                for j=1:npoints

                    xu = m{cam}(1,j);
                    yu = m{cam}(2,j);
                    r = sqrt(xu^2 + yu^2);
                    
                    rpoly = polyval(coefs, r);
                    
                    m{cam}(1,j) = xu * rpoly;
                    m{cam}(2,j) = yu * rpoly;
                end
            end
            
            
        elseif (modeltype > 0)

            for cam=1:ncams

                if size(camcalib, 2) > 1
                    
                    K = camcalib(((i-1)*3+1):(i*3) ,1:3);
                    Ki = inv(K);
                else
                    K = eye(3);
                    Ki = eye(3);
                end

                
                if size(rdcoefs, 2) == ncams

                    lambda = rdcoefs(cam);
                else
                    lambda = rdcoefs(1);
                end       

                if lambda == 0 
                    continue;
                end

                for j=1:npoints

                    xu = m{cam}(1,j);
                    yu = m{cam}(2,j);
                    
                    xM = m3d{cam}(:,j);
                    xM = Ki*xM;
                    nm = norm(xM);
                    
                    fi = atan2(xM(1), xM(2));
                    theta = acos(xM(3) / nm);
                        
                    switch (modeltype)
                        
                        case 1 % stereographic projection r = k tan(theta)
                            r = lambda * tan(theta);
                            
                        case 2 % equidistant projection r = ktheta,
                            r = lambda * theta;

                        case 3 % equisolid angle projection r = k sin theta/2
                            r = lambda * sin(theta/2);
                        
                        case 4 % sine law projection r = k sin theta
                            r = lambda * sin(theta);
                    end
                    
                    p = K*[r * cos(fi); r * sin(fi); 1];
                    p = p/p(3);
                                                             
                    m{cam}(1,j) = p(1);
                    m{cam}(2,j) = p(2);
                end
            end
            
        end
    end

    %add outliers
    noutlierscnt = floor(noutliers * npoints);
    perm = randperm(npoints);
    outs=ones(3, noutlierscnt);
    for j=1:ncams

        if size(camcalib, 2) > 1
            K = camcalib(((j-1)*3+1):(j*3), 1:3);
        else
            K = eye(3);
        end
        %outs = radius(1)*(rand(3, noutlierscnt) - 0.5);
        outs(1,:) = rand(1, noutlierscnt)*(screenSize(1,2)-screenSize(1,1)) + screenSize(1,1);
        outs(2,:) = rand(1, noutlierscnt)*(screenSize(2,2)-screenSize(2,1)) + screenSize(2,1);
%        outsimg = K*outs;
%        outsimg = outsimg ./ repmat(outsimg(3,:), 3, 1);
        m{j}(:,perm(1:noutlierscnt))=outs(1:2,:);
    end
    
    if length(anoise) < 2
        anoise = anoise*ones(1,ncams);
    end
    
    %add 2D noise
    for j=1:ncams

        % p = m{j} + (rand(2,npoints) - 0.5)*anoise;
        p = m{j} + (1/3*randn(2,npoints))*anoise(j);
        m{j}=p;
    end
    
end

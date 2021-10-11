% show camera and 3D content
%
% by Martin Bujnak

% Pgt - 3x4 matrices in cell array
% Kgt - 3x3 calibration matrices in cell array (or single K for all equal)
% m - cells with 3xN arra of 2D points
% Mgt - 3xN array of 3D points
%
% show2D - show figures with reprojections of 3D points to 2D
%   showMotion - display arrows from corresponding points in previous
%   camera
% show3D - show cameras in 3D

function [fig_number] = ShowCameras(Pgt, Kgt, m, Mgt, show2D, showMotion, show3D, markpoints, m2points)
   
    fig_number = [];
    camCnt = length(m);
    
    if ~iscell(Kgt)
        Kgt={Kgt};
    end

    % show cameras...
    if show2D
        for i=1:camCnt

            % display 2D features
            fig_number = figure;

            set(fig_number,'Units','normalized', ...
                'BackingStore','off', ...
                'Color',[0.8 0.8 0.8], ...
                ...%'MenuBar','none', ...
                'Resize','on', ...
                'Name', 'flow', ...
                'Position',[0 1-1/2 1/3 1/2-1/25], ...
                'NumberTitle','off'); 

            xlabel('X-axis');
            ylabel('Y-axis');
            xlim([-1,1]);
            ylim([-1,1]);
            box; 
            set(gca,'dataaspectratio',[1 1 1])
            hold on;

            plot(m{i}(1,:), m{i}(2,:), '.r');
            
            if nargin > 8
                plot(m2points{i}(1,:), m2points{i}(2,:), '.g');
            end

            if showMotion

                if i == camCnt
                    j = 1;
                else 
                    j = i + 1;
                end

                plot(m{j}(1,:), m{j}(2,:), '.g');
                quiver(m{j}(1,:), m{j}(2,:), m{i}(1,:)-m{j}(1,:), m{i}(2,:)-m{j}(2,:), 0, 'g--');
            end        
            
            plot(m{i}(1,markpoints), m{i}(2,markpoints), 'ob');
        end
    end    
    
    if show3D
        
        % show 3D scene
        figure;
        axis equal
        hold on;
        plot3(Mgt(1,:), Mgt(2,:), Mgt(3,:), '.g');

        plot3(Mgt(1,markpoints), Mgt(2,markpoints), Mgt(3,markpoints), 'ob');
        
        calibCnt = length(Kgt);

        % show cameras...
        for i=1:camCnt

            if calibCnt > 1

                P = inv(Kgt{i})*Pgt{i};
            else

                P = inv(Kgt{1})*Pgt{i};
            end

            R = inv(P(1:3,1:3));
            pos = -R*P(:,4);

            plot3(pos(1), pos(2), pos(3), 'o');

            R1 = 15*R(:,1);
            R2 = 15*R(:,2);
            R3 = 15*R(:,3);

            plot3(pos(1)+[0 R1(1)], pos(2)+[0 R1(2)], pos(3)+[0 R1(3)], '-r');
            plot3(pos(1)+[0 R2(1)], pos(2)+[0 R2(2)], pos(3)+[0 R2(3)], '-g');
            plot3(pos(1)+[0 R3(1)], pos(2)+[0 R3(2)], pos(3)+[0 R3(3)], '-b');

        end
    end
    
end
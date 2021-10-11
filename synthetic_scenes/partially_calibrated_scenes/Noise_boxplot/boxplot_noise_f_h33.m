name = '';


%number of methods
nmet = 2;

%noise levels
fnoises = noise_levels;

% methods names
methodsNames = { 'sH5f_{2}' 'sH5f_{3}' };
% name = 'Forward motion';

% methods for displaying
methods = [1 2];

%limits of y-axis (can be read from the data)
ymin = -14; ymax = -12;

% colors for methods
cols = {'r' 'b'};


errors_to_plot = {};
for j = 1:nmet
    for ns=1:length(fnoises)   
        errors_to_plot{j}(ns).err = (errors{j}(ns).f_err);
    end
end


show_boxplot(errors_to_plot, fnoises, methodsNames, name, methods, ymin, ymax, cols)



function show_boxplot(algs, fnoises, methodsNames, name, methods, ymin, ymax, cols)


ncnt = length(fnoises);

%  allmes{cfg}{focal}{alg}(noise, measurement)
i = 1;
ii = 1;


hnd = figure;
axes_handle=axes('fontsize', 24);
hold on;
set(hnd, 'Name', name);
title(name, 'fontsize', 20);
lims =  [ymin   ymax];
ylim(axes_handle, lims);
xlim(axes_handle, [0.4 length(fnoises)+0.5]);

xt = [];
xn = {};
for ns=1:length(fnoises)
    xt(ns) = ns*1;
    xn{ns} = num2str(fnoises(ns));
    
end

for ns=1:length(methods)
    plot(-1,-1, cols{ns});
end

set(gca,'XTick', xt);
set(gca,'XTickLabel', xn);

xlabel('Noise \sigma (pixels)', 'fontsize', 24);
% xlabel('Distance of 3D point from plane', 'fontsize', 24);
% ylabel('$Log_{10} (|| t - t_{gt} ||_{2}/|| t_{gt} ||_{2})$','Interpreter','Latex');
ylabel('Relative focal length error(%)','fontsize', 24);

for ns=2:length(fnoises)
    plot([ns-0.5 ns-0.5], [lims(1) lims(2)], ['-'], 'LineWidth', 2, 'color', [0.9 0.9 0.9]);
end

metshift = 0.9/length(methods);
xsize = metshift/2-0.02;
xxsize = metshift/4-0.02;
ofs = (1 + xsize) - (metshift * length(methods)) / 2;



%plot([0 10], [1 1], 'Color', 'cyan');

scollector = [];
for j=methods
    
    xstart = ofs;
    xpos = xstart;
    
    for ns=1:length(fnoises)
        
        x = [];
        
        % error that will be plot
        %   mx =abs(algs{j}(ns).f(:,1)-(f1gt*pixel))/(f1gt*pixel);
        mx =algs{j}(ns).err;
        x = [x; mx(:)];
        
        
        
        [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = mb_boxutil(x, 1,1.5,0);
        
        plot([xpos xpos], [q1 n2], ['--' cols{ii}], 'LineWidth', 2);
        plot([xpos xpos], [q3 n1], ['--' cols{ii}], 'LineWidth', 2);
        
%         plot(xpos * ones(1, length(yy)), yy, ['+' cols{ii}], 'MarkerSize', 3);
        
        plot([xpos-xxsize xpos+xxsize], [q1 q1], ['-' cols{ii}], 'LineWidth', 2);
        plot([xpos-xxsize xpos+xxsize], [q3 q3], ['-' cols{ii}], 'LineWidth', 2);
        
        plot([xpos-xsize xpos+xsize xpos+xsize xpos-xsize xpos-xsize], [n1 n1 n2 n2 n1], ['-' cols{ii}], 'LineWidth', 2);
        
        plot([xpos-xsize xpos+xsize], [med med], ['-' cols{ii}], 'LineWidth', 2);
        
        stat = [n2, med, n1];
        
        %scollector(:, i) = stat;
        %collector(:, i) = x;
        i=i+1;
        xpos = xpos + 1;
    end
    
    %X = (1:size(lambdas,2)) + xstart -1;
    
    
    ofs = ofs + metshift;
    ii = ii + 1;
 
%     legend(methodsNames(methods), 'Location','northwest');
    
end
end

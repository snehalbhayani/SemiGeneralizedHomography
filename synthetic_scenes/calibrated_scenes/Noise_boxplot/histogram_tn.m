%number of methods
nmet = 8;

% methods names
methodsNames = { 'sH5_{2}' 'sH4.5_{2}'  'sH5_{3}' 'sH4.5_{3}'  'sH5_{4}'  'P3P+N', 'E_{5+1}', 'E_{4+2}'};
name = 'Translation error';


% methods for displaying
methods = [1 2 3 4 5 6 7 8];

cfgstyles = {'-' '-.*' '--' '--' ':'};

%limits of y-axis (can be read from the data)
xmin = -20;
xmax = 5;

% colors for methods
cols = {'r' ,[0.9648    0.8086    0.4141],'b' ,[0.3906    0.6992    0.6055], 'k' ,  'g',  'c', 'm' };

errors_to_plot = {};
for j = 1:nmet
    
    errors_to_plot{j}.err = (errors{j}.tn_err);
    
end




show_boxplot(errors_to_plot, methodsNames, name, methods, xmin, xmax, cols, cfgstyles)



function show_boxplot(errors_to_plot, methodsNames, name, methods, xmin, xmax, cols, cfgstyles)


hnd = figure;
axes('fontsize', 24);
hold on;
set(hnd, 'Name', name);

%title(['focal length 1 ' int2str(fc*18) 'mm']);
% title(name, 'fontsize', 12);
% noisefree

% xlabel('$Log_{10} (|| t - t_{gt} ||_{2}/|| t_{gt} ||_{2})$','Interpreter','Latex');
xlabel('Log_{10} relative pose error ', 'fontsize', 24);
ylabel('Frequency', 'fontsize', 24);

log_range = xmin:0.3:xmax;
cfgs = [1];
% log10 bez 0
for cfg=cfgs
    
    ii = 1;
    for j=methods
        
        data = [];
        
        data = [data  errors_to_plot{j}.err];
        
        
        if 1
            gr = data;
            
            % filter exact data
            used = find(gr ~= 0);
            totalmin = 1e-20;
            used = find(gr == 0);
            gr(used) = totalmin;
            
            h = hist(log10(abs(gr)), log_range);
            
            % smoothing
            if 1
                ks = 3;
                smask = [0.2 0.5 1 2 1 0.5 0.2];
                %s = 2;
                %smask = [fliplr(exp(-(1:ks).^s)) 1 exp(-(1:ks).^s)];
                v = conv(h, smask)*(1/sum(smask));
                h = v(ks:end-(ks+1));
            end
            
            plot( log_range , h, [cfgstyles{cfg}],'color', cols{ii}, 'LineWidth', 2);
            ii = ii + 1;
        end
    end
end

% add legend
legend(methodsNames(methods), 'fontsize', 20, 'Orientation','horizontal', 'NumColumns', 2, 'NumColumnsMode', 'manual');
end


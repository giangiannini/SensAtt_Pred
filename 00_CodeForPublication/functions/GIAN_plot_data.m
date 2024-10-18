function handle_figure = gian_plot_data(cfg, varargin)
    %this function created by Gian is used to plot the ERPs in a similar
    %way fieldtrip does but with Variance bars and 
    for i = 1:length(varargin)
        y{i} = varargin{i}.avg;
        x{i} = varargin{i}.time;
        err{i} = varargin{i}.err;
        if isfield(cfg, 'baseline')
            samples_baseline = 2:24;
            y{i} = y{i} - mean(y{i}(1:64,samples_baseline), [2]);
        end
    end

    %set(gca, 'XScale', 'log')
    chan = []; 
    if isstring(cfg.channel)
        for i = 1:length(cfg.channel)
            chan = [chan find(strcmp(varargin{1}.label, cfg.channel(i)))];
        end
    else
        chan = cfg.channel; 
    end

    %% DO LAYERS
    if ~isfield(cfg, 'printlayers')
        cfg.printlayers = 0; 
    end

    if cfg.printlayers == 1
        for i = 1:length(varargin)
            figure('Renderer', 'painters', 'Position', [10 10 900 600])
            l(i) = plot(x{i}, mean(y{i}(chan,:),1), 'Color', cfg.colors(i,:), 'LineStyle', cfg.linestyle{i});
            hold on
            % if length(cfg.colors(i,:)) > 3
            %     set(l, 'FaceAlpha', cfg.colors(i,4))
            % end
            if isfield(cfg, 'xlim')
                xlim(cfg.xlim)
            end
            if isfield(cfg, 'ylim')
                ylim(cfg.ylim)
            end

            p = patch([x{i} flip(x{i})], [mean(y{i}(chan,:),1)-mean(err{i}(chan,:),1) flip(mean(y{i}(chan,:),1)+mean(err{i}(chan,:),1))], ...
                'b', 'FaceAlpha',0.25, 'EdgeColor','none');
            set(p, 'FaceColor', cfg.colors(i,[1:3]))

            %plot title
            if length(cfg.channel) == 1
                z = ylim;
                txt = char(cfg.channel);
                t = text(0.40,z(2)-0.45,txt);
                t.FontSize = 25;
            elseif length(cfg.channel) > 1
                cfg.channel = string(cfg.channel); 
                plot_title = strcat("mean(", strjoin(cfg.channel, ', '), ")");
                title(plot_title)
            end
            exportgraphics(gcf, strcat(cfg.output_printlayers, '_', num2str(i), '_', '.emf'))
            close all
        end

        %PLOT ALSO SIGNIFICANCE AREA
        if isfield(cfg, 'mask')
            figure('Renderer', 'painters', 'Position', [10 10 900 600])
            if isfield(cfg, 'xlim')
                xlim(cfg.xlim)
            end
            if isfield(cfg, 'ylim')
                ylim(cfg.ylim)
            end

            [y_pos x_pos] = find(cfg.mask(chan,:)==1);
            bands = [x{1}(min(x_pos)) x{1}(max(x_pos))]; 
            hold on
            xp = [bands fliplr(bands)];                                                         % X-Coordinate Band Definitions 
            yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';                  % Y-Coordinate Band Definitions
            for k = 1:size(bands,1)                                                             % Plot Bands
                patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.25, 'EdgeColor', 'none')
            end

            %plot title
            if length(cfg.channel) == 1
                z = ylim;
                txt = char(cfg.channel);
                t = text(0.40,z(2)-0.45,txt);
                t.FontSize = 25;
            elseif length(cfg.channel) > 1
                cfg.channel = string(cfg.channel); 
                plot_title = strcat("mean(", strjoin(cfg.channel, ', '), ")");
                title(plot_title)
            end
            exportgraphics(gcf, strcat(cfg.output_printlayers, '_significance_', '.emf'))
            close all
        end
    end


    %% SHOW PROPER IMAGE
    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    for i = 1:length(varargin)
        l(i) = plot(x{i}, mean(y{i}(chan,:),1), 'Color', cfg.colors(i,:), 'LineStyle', cfg.linestyle{i});
        hold on
        % if length(cfg.colors(i,:)) > 3
        %     set(l, 'FaceAlpha', cfg.colors(i,4))
        % end
    end
    if isfield(cfg, 'xlim')
        xlim(cfg.xlim)
    end
    if isfield(cfg, 'ylim')
        ylim(cfg.ylim)
    end

    for i = 1:length(varargin)
        p = patch([x{i} flip(x{i})], [mean(y{i}(chan,:),1)-mean(err{i}(chan,:),1) flip(mean(y{i}(chan,:),1)+mean(err{i}(chan,:),1))], ...
               'b', 'FaceAlpha',0.25, 'EdgeColor','none');
        set(p, 'FaceColor', cfg.colors(i,[1:3]))
    end
    hold off

    if isfield(cfg, 'legend')
        legend([l(1:i)], cfg.legend(1:i))
    end


    %plot mask
    if isfield(cfg, 'mask')
        [y_pos x_pos] = find(cfg.mask(chan,:)==1);
        bands = [x{1}(min(x_pos)) x{1}(max(x_pos))]; 
        hold on
        xp = [bands fliplr(bands)];                                                         % X-Coordinate Band Definitions 
        yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';                  % Y-Coordinate Band Definitions
        for k = 1:size(bands,1)                                                             % Plot Bands
            patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.25, 'EdgeColor', 'none')
        end
    end

    %plot title
    if length(cfg.channel) == 1
        z = ylim;
        txt = char(cfg.channel);
        t = text(0.40,z(2)-0.45,txt);
        t.FontSize = 25;
    elseif length(cfg.channel) > 1
        cfg.channel = string(cfg.channel); 
        plot_title = strcat("mean(", strjoin(cfg.channel, ', '), ")");
        title(plot_title)
    end
        
    hold off
    xlabel('s')
    ylabel('mV')

end
function whole_time_plot(t, x, p, colours)

fh = figure();
clf()
for i = 0:p.num_lanes-1
    data = x(:,(i*p.len_lane+3):((i+1)*p.len_lane-2)).';
    subplot(p.num_lanes,1,i+1); 
    hm = heatmap(t,1:p.cell_per_lane,data);
    yname = sprintf('Lane %d',i);
    if i ~= p.num_lanes-1
        Ax = gca;
        Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    end
    ylabel(yname);
    colormap(colours(1:p.number_species,:)); %number of species
    
    % Why is the plotting here incorrect though? 
    
    % Put an axis over top of the heatmap
    % The axis will not be visible but will cover all area 
    % to the right of the heatmap. 
    hmp = hm.Position; 
    cbax = axes('Position',[sum(hmp([1,3])), 0, 1-sum(hmp([1,3])), 1],...
        'XTick',[], 'YTick', [], 'Color',fh.Color);
    cbax.XAxis.Visible = 'off';
    cbax.YAxis.Visible = 'off';

    % Set the new axis color map to match the 
    % heatmap's colormap
    cbax.Colormap = hm.Colormap; 

    % Add a colorbar the same vertical position as the heatmap
    cbh = colorbar(cbax,'Position',[.90, hmp(2), .03, hmp(4)]); 
    % Set the limits to 0:1 and set the ticks so that they 
    % are in the center of each bar
    cbh.Limits = [0,1]; 
    nColors = size(hm.Colormap,1); 
    cbh.Ticks = (1/nColors : 1/nColors : 1) - 1/nColors/2;

    % Set the new labels for each tick
    cbh.TickLabels = 0:nColors-1;
end
xlabel('time (s)');
saveas(gcf,'lane_history.pdf')

end

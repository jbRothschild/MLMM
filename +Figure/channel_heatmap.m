function fh = channel_heatmap(time_point, x, p, colours, show)

fh = figure()
if ~show
   set(fh,'Visible','off');
   
clf()
x0=10;
y0=10;
width=p.cell_per_lane*125;
height=p.num_lanes*50;
set(gcf,'position',[x0,y0,width,height])

x_at_t = zeros(p.num_lanes,p.cell_per_lane);
for i = 0:p.num_lanes-1
    x_at_t(i+1,:) = x(time_point,(i*p.len_lane+3):((i+1)*p.len_lane-2));
end
hm = heatmap(1:p.cell_per_lane,1:p.num_lanes,x_at_t);
colormap(colours(1:p.number_species,:)); %number of species
ylabel('Lane');
xlabel('Cell position');

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
%% Which simulation to run
%sfun = @simulation
sfun = @simulation_nbr_cells_per_lane;
nbr_simulations = 5000;

%% Default parameters
parameters = cell(3,1);

init       = zeros(2,5);
init(1,:)  = ones(1,size(init,2));

%                  Green  Red   ???
parameters{1,1} = [ 0.01  0.01   % up-left
                   0.1   0.1    % left
                   0.01  0.01   % down-left
                   0.01  0.01   % up-right
                   0.01  0.01   % right
                   0.01  0.01   % down-right
                  ];      % reactions
parameters{2,1} = [0, 1000000]; % time (seconds)
parameters{3,1} = init; % init   % 


%% Function Call
sfun(parameters, nbr_simulations);

%% Nothing changes
function [mean_fpt] = simulation(parameters, nbr_sims, varargin)
import Figure.*
nbr_species = size(parameters{1,1},2);
fpt    = zeros(1,nbr_sims);
winner = zeros(1,nbr_sims);
init   = parameters{3,1};
final  = zeros(size(init,1),size(init,2));
    
for i=1:nbr_sims
    [t,x,p,colours]  = GillespieMLMM(parameters{1,1},parameters{2,1},parameters{3,1});
    fpt(1,i) = t(end);
    for j=0:size(init,1)-1
        final(j+1,:) = x(end,j*(size(init,2)+4)+3: j*(size(init,2)+4)+2+size(init,2));
    end
    winner(1,i) = max(final,[],'all');   
end

if isempty(varargin)
    figure(1)
    %{
    s1 = subplot(2,2,1);
    hold on
    fig1 = channel_heatmap(1, x, p, colours, false);
    copyobj(s1,fig1);
    hold off

    s2 = subplot(2,2,2);
    hold on
    for i=0:nbr_species-1
        sum(winner(:) == i)/nbr_sims
        bar([i],[sum(winner(:) == i)/nbr_sims],'FaceColor',colours(i+1,:),1)
    end
    hold off

    s3 = subplot(2,2,[3,4]);
    hold on
    %}
    [N,edges] = histcounts(fpt);
    centers   = (edges(1:end-1) + edges(2:end))./2;
    scatter(centers,N,10,'filled')
    xlabel('count');
    ylabel('time (s)');
    title('fpt of takeover');
    
    saveas(gcf,'fpt_default_sim.pdf');
end

mean_fpt = mean(fpt);

end

%% Varying number of lanes
function simulation_nbr_lanes(parameters, nbr_sims)

nbr_cells_per_lane = 10;
nbr_lanes          = 10;

mean_fpt = zeros(1,nbr_lanes);

for i=2:nbr_cells_per_lane
    % TODO
    
    mean_fpt(i-1) = simulation(parameters, nbr_sims, false); 
    
end  
scatter(2:nbr_cells_per_lane,mean_fpt,10,'filled')
xlabel('Number of lanes');
ylabel('Mean fpt (s)');
title('Varying number lanes');

end

%% Varying number of cells in each lane
function simulation_nbr_cells_per_lane(parameters, nbr_sims)

nbr_cells_per_lane = 20;
nbr_lanes          = 2;

mean_fpt = zeros(1,nbr_cells_per_lane-2);

for i=3:nbr_cells_per_lane
    init       = zeros(nbr_lanes,i);
    init(1,:)  = ones(1,size(init,2));
    parameters{3,1} = init;
    
    mean_fpt(i-2) = simulation(parameters, nbr_sims, false);
    
end  

figure()
scatter(3:nbr_cells_per_lane,mean_fpt,10,'filled')
xlabel('cells per lane');
ylabel('mean fpt (s)');
title('Varying lane length');

saveas(gcf,'fpt_vary_lane_length.pdf');

end

%% Number of species change
function simularion_nbr_species(paramaters, nbr_sims)

end
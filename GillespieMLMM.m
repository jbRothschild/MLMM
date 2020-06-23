function [t,x,p,colours] = GillespieMLMM(varargin) % TODO : Get rid of p, colours
% MMLM simulation with simple Gillespie
import Gillespie.*
%% Reaction network:
%   For now, each of the cells has 6 reactions, it can push out a neighbour
%   in each direction and switch to another lane, puching the neighbor out.
%  Rate constants
%   It's possible that the different types of cells have their own rates,
%   in which case you would define them here. FOr now they're the same,
%   with [Green Red], or [0 1]
%   The scheme below defines, for each cell, what the reacion is they do.
%   Note that the column is the cell type, row is the cell reaction
%   direction

if isempty(varargin)
    % Default
    %             Green  Red 
    p.cell_rxn = [ 0.01  0.01   % up-left
                   0.1   0.1    % left
                   0.01  0.01   % down-left
                   0.01  0.01   % up-right
                   0.01  0.01   % right
                   0.01  0.01   % down-right
                  ];      
    tspan = [0, 1000000]; % seconds
    init  = randi([0 size(p.cell_rxn,2)-1], 2, 5)     % 0 -> Green      1 -> Red        i -> Next species

else
    % Import values
    p.cell_rxn = cell2mat(varargin(1));
    tspan      = cell2mat(varargin(2));
    init       = cell2mat(varargin(3));
end

p.num_lanes      = size(init,1);
p.cell_per_lane  = size(init,2);
p.number_species = size(p.cell_rxn,2);
p.len_lane       = p.cell_per_lane + 4; % empty cell at each end
p.tot_num_cells  = p.num_lanes*p.cell_per_lane;
p.tot_positions  = p.num_lanes*p.len_lane;

% transform init to proper form
x0               = zeros(1,p.tot_positions);
for i=0:p.num_lanes-1
    x0(i*p.len_lane+3: i*p.len_lane+2+p.cell_per_lane) = init(i+1,:);
end
%% Specify reaction network
pfun = @propensities_MLMM;    % defines your probabilities
sfun = @stoichiometries_MLMM; % defines your updates

%% Run simulation
[t,x] = directMethodMLMM(sfun, pfun, tspan, x0, p);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);

%% Plot time course
% TODO : It's going to be a bit of a challenge to figure out what plot to
%        show. We have a couple of lanes, so it might seem more useful to
%        plot each lane in time? Or a .gif/animation of the whole thing?
colours = [34./255., 139./255., 34./255. ; 139./255., 0, 0 ; 0, 0, 139./255. ; 0, 139./255., 0. ];
if isempty(varargin)
    import Figure.*
    colours = [34./255., 139./255., 34./255. ; 139./255., 0, 0 ; 0, 0, 139./255. ; 0, 139./255., 0. ];
    channel_heatmap(1, x, p, colours, false)
    whole_time_plot(t, x, p, colours)
    gif_channel(t, x, p, colours)
end


end

%% Stoichiometries
function x_next = stoichiometries_MLMM(x, mu, p)
% Return the next state of the system
%    x   :   current state of the system
%    mu  :   reaction number
%    p   :   parameters
%    
%    -To get the lane number: floor(cell_number/len_lane) + 1
%    -To get the position of the cell in the lane:
%     mod(cell_number,len_lane)

rxn_per_cell = size(p.cell_rxn,1);
total_rxn = p.tot_positions * rxn_per_cell;

cell_number = floor((mu-1)/rxn_per_cell) + 1;
cell_rxn_number = mod((mu-1),rxn_per_cell) + 1;
cell_lane = floor(cell_number/p.len_lane) + 1; % lane of the cell
cell_position_in_lane = mod(cell_number,p.len_lane);

if ismember( cell_rxn_number,[1 2 3] ) % Moving left
    lane_change = cell_rxn_number - 2; % Whether the lane is changed -1,0,1
    cell_to_move = cell_number + lane_change*p.len_lane - 1; % cell that moves the rest
    cell_to_move_position_in_lane = cell_position_in_lane - 1; % position in lane
    x( (cell_to_move - cell_to_move_position_in_lane+1) : (cell_to_move-1) ) ...
            = x( (cell_to_move - cell_to_move_position_in_lane+2):(cell_to_move) );
    
elseif ismember( cell_rxn_number,[4 5 6] ) % Moving right
    lane_change = cell_rxn_number - 5; % Whether the lane is changed -1,0,1
    cell_to_move = cell_number + lane_change*p.len_lane + 1; % cell that moves the rest
    cell_to_move_position_in_lane = cell_position_in_lane + 1;
    x( (cell_to_move+1):(cell_to_move + p.len_lane - cell_to_move_position_in_lane ) ) ...
            = x( (cell_to_move):(cell_to_move + p.len_lane - cell_to_move_position_in_lane - 1) );

% Cell that split
x(cell_to_move) = x(cell_number);
        
    
end

x_next = x;

end

%% Propensities
function a = propensities_MLMM(x, p)
% Return reaction propensities given current state x
%    -In these reactions, to get the cell number of the reaction, simply
%     divide floor(tot_rxn/rxn_per_cell) + 1
%    -To get the lane number: floor(cell_number/len_lane) + 1
%    -To get the position of the cell in the lane:
%     mod(cell_number,len_lane)
rxn_per_cell = size(p.cell_rxn,1);
total_rxn = p.tot_positions * rxn_per_cell;

a = zeros(1,total_rxn);
for i = 1:size(a,2)
    % cell_number : floor((i-1)/rxn_per_cell) + 1
    % rxn_number  : mod((i-1),rxn_per_cell) + 1
    % cell_type : x(cell_number) + 1. i.e. [0 1] -> [1 2] for labeling 
    a(i) = p.cell_rxn(mod((i-1),rxn_per_cell) + 1, ...           
                      x( floor((i-1)/rxn_per_cell) + 1 ) + 1 );
end

% if in side lane, can't move into the wall, so set switching to 0
for i = 0:(p.len_lane-1)
    a(i*rxn_per_cell + 1) = 0.0;
    a(i*rxn_per_cell + 4) = 0.0;
    a(p.len_lane*(p.num_lanes-1)*rxn_per_cell + i*rxn_per_cell + 3) = 0.0;
    a(p.len_lane*(p.num_lanes-1)*rxn_per_cell + i*rxn_per_cell + 6) = 0.0;
    
end

% positions outside of the channel are dead
for i = 0:(p.num_lanes-1)
    for j = 1:rxn_per_cell
        a( i*p.len_lane*rxn_per_cell + j ) = 0.0;
        a( i*p.len_lane*rxn_per_cell + rxn_per_cell + j ) = 0.0;
        a( (i+1)*(p.len_lane)*rxn_per_cell + 1 - j ) = 0.0;
        a( (i+1)*(p.len_lane)*rxn_per_cell + 1 - rxn_per_cell - j ) = 0.0;
    end
end
end





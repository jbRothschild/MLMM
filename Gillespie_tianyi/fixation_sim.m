row = 4;
column = 10;
grid = zeros(row, column);   %create a grid to put cells 
grid(1,:) = 1;               %"1" represents a green cell, "0" represents a red cell
grid(2,:) = 1;

% All the rates are in minute^-1

%rates for red cells
r1 = 0.1;            %rate of the cell stay in the lane
r2 = 0.01;           %rate of lane changing
r = r1+r2;           %Sum of all rates 
r_s = r1 /r;         %probability of red cell stay in the lane 
                     %probability of red cell change lane is 1-r_s

%rates for gree cells
g1 = 0.1;
g2 = 0.01;
g = g1+g2;
g_s = g1 /g;


t = 0;               %initial time
T = 10;             %total time 
Cells_0 = grid;        %inital cell grids 
Cellgird = cell(numel(T),1);
Cellgrid{1} = Cells_0;
count = 1;
while t<T
    Cells = [Cellgrid{count}];
    if all(Cells>0.5) | all(Cells<0.5)   %run the loop till one type of cells dominate the entire grid 
        break
    end
    
    tau_r = -log(rand(size(Cells)))/r;     %random numbers from exponential distribition to entire matix 
    tau_r = tau_r + Cells;                 %since red cells are zeors and green are ones, this can return the minmum of red cells
    [t_r, mu_r] = min(tau_r,[],2);         %mu_r gives all the colume number of the minimum in each row 
    [dt_r,I_r] = min(t_r);                 %I_r gives the row number of the minimum of entir matix
                                           %dt_r gives the minimum time 
            
    tau_g = -log(rand(size(Cells)))/g;    
    tau_g = tau_g - Cells;                  
    [t_g, mu_g] = min(tau_g,[],2);
    [dt_g, I_g] = min(t_g);                %dt_g gives the minimum time of green cells 
    
    if dt_r < (dt_g+1)                     %Compare green and red cells' minimum time and pick the smaller one 
        dt = dt_r;
        I = [I_r, mu_r(I_r)];              %Locate the cell has the minimum time    
    else 
        dt = dt_g+1;
        I = [I_g, mu_g(I_g)];
    end
    
    celltype = Cells(I(1),I(2));           %Find the colour of the cell that has the minimum time: 1 -- green; 0 -- red
    d = 4;                                 %number of lane changing directions (diagonal)
    h = 2;                                 %number of growth directions (horizontal)
    
    if celltype == 1                       %check the colour to use the relative rates 
        stay = g_s;
    else 
        stay = r_s;
    end
    
    [new_row, new_column] = Direction(stay, d, h, I);
    
    if new_row > 0 && new_row < row         %aviod the problem from cells at the edges
        if new_row == I(1)                   %stay in lane 
            if new_column < I(2)             %shift to the left
                if new_column > 0
                    for i = 1: new_column
                        Cells(I(1), i) = Cells(I(1), i+1);
                    end
                end
            else                            %shift to the right
                if new_column < column  
                     for i = column:-1:new_column
                        Cells(I(1), i) = Cells(I(1), i-1);
                    end
                end
            end
        else 
            if new_column <I(2)             %top/bottom-left
                if new_column > 0
                    for i = 1:new_column
                        Cells(new_row,i)= Cells(new-row, i+1);
                    end
                    for i = I(2): column
                        Cells(I(1),i) = Cells(I(1),i+1);
                    end
                end
            else                            %top/bottom-right
                if new_column < column
                    for i = column:-1:new_column
                        Cells(new_row,i)= Cells(new-row, i-1);
                    end
                    for i = I(2):-1:1
                        Cells(I(1),i) = Cells(I(1),i-1);
                    end
                end
            end
        end 
    end 
    
    t = t + dt;
    count = count+1;
    Cellgrid{count} = Cells;
end

function [new_row, new_column] = Direction(stay,d,h, I)
    if rand()< stay                         %check if stay in the lane or not
        if rand() < 1/h                %check the dirction 
            new_row = I(1);            %find the location of next state
            new_column = I(2)-1;
        else 
            new_row = I(1);  
            new_column = I(2)+1;
        end
     else
        if rand() <1/d
            new_row = I(1)-1;
            new_column = I(2)-1;                 
        elseif rand()<2/d && rand()>=1/d  
            new_row = I(1)-1;
            new_column = I(2)+1;                     
        elseif rand()<3/d && rand()>=2/d  
            new_row = I(1)+1;
            new_column = I(2)-1;                  
        else 
            new_row = I(1)+1;
            new_column = I(2)+1;
        end
    end
end


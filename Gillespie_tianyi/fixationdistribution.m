clear all;
close all;

runs = 100;
Winner = zeros(1,runs);
Fixationtime = zeros(1,runs);
for n = 1:runs
    [winner , time] = Fixation(4,10,[3,4], [1,2], 1000, [0.1,0.01,0.1,0.01]);
    Winner(n) = winner;
    Fixationtime(n) = time;
end

notfix = find(Winner ==3);
Winner(notfix)=[];
Fixationtime(notfix)=[];

red = find(Winner == 0);
redtime = Fixationtime(red);

green = find(Winner == 1);
greentime = Fixationtime(green);

figure
pie([sum(red),sum(green),runs])
labels = {'Red wins','Green wins','No winner'};
legend(labels,'Location','southoutside','Orientation','horizontal','FontSize',20);

hg = histogram(greentime);
%hg.Normalization = 'probability';
hg.BinWidth= 25;
hold on;
hr = histogram(redtime);
%hr.Normalization = 'probability';
hr.BinWidth = 25;
xlabel('Fixation time (minutes)','FontSize', 20);
ylabel('Counts','FontSize', 20);
title('500 runs: r1 = 0.1/min, g1 = 0.1/min (growth)   r2 = 0.01/min, g2 = 0.01/min (lane changing)','FontSize', 20);
legend('Green cells', 'Red cells', 'FontSize', 15);
                 


function [winner, time] = Fixation(row,column,redlane, greenlane,totaltime,rates)
   
grid = zeros(row, column);
for i = redlane
    grid(i,:) = 0;
end
for j = greenlane
    grid(j,:)= 1;
end
Cellgrid{1} = grid;

T = totaltime;
t = 0;
count = 1;

r = rates(1)+rates(2);
r_s = rates(1)/r;

g = rates(3)+rates(4);
g_s = rates(3)/g;

while t<T
    Cells = [Cellgrid{count}];
    check = Cells(:);
    if sum(check == 1) == (row*column) | sum(check==0) == (row*column) %run the loop till one type of cells dominate the entire grid 
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
    
    [new_row, new_column] = Direction(stay, row, d, h, I);
    
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
        oldcell = Cells(I(1), I(2));
        if new_column <I(2) && new_column > 1             %top/bottom-left
            
            for i = 1:new_column-1
                Cells(new_row,i)= Cells(new_row, i+1);
            end
            Cells(new_row, new_column) = oldcell;
        elseif new_column >I(2) && new_column <column          %top/bottom-right
            
            for i = column:-1:new_column+1
                Cells(new_row,i)= Cells(new_row, i-1);
            end
            Cells(new_row, new_column) = oldcell;
        end
        
    end 
    
    t = t + dt;
    count = count+1;
    Cellgrid{count} = Cells;
    timestep(count) = t;
end

final = Cellgrid{end};
if all(final(:) == 1)
    winner = 1; 
elseif all(final(:) == 0)
    winner = 0;
else
    winner = 3;
end

time = t;
end

function [new_row, new_column] = Direction(stay,row,d,h, I)
    if rand()< stay                         %check if stay in the lane or not
        if rand() < 1/h                %check the dirction 
            new_row = I(1);            %find the location of next state
            new_column = I(2)-1;
        else 
            new_row = I(1);  
            new_column = I(2)+1;
        end
    else
         if I(1) == 1
             d = 2;
             if rand()< (1/d)
                 new_row = I(1)+1;
                 new_column = I(2)-1;
             else 
                 new_row = I(1)+1;
                 new_column = I(2)+1;
             end
         elseif I(1) == row
             d = 2;
             if rand()< (1/d)
                 new_row = I(1)-1;
                 new_column = I(2)-1;
             else 
                 new_row = I(1)-1;
                 new_column = I(2)+1;
             end
         else 
             if rand() <(1/d)
                 new_row = I(1)-1;
                 new_column = I(2)-1;                 
             elseif rand()< (2/d) && rand()>=(1/d) 
                 new_row = I(1)-1;
                 new_column = I(2)+1;                     
             elseif rand()<(3/d) && rand()>=(2/d)  
                 new_row = I(1)+1;
                 new_column = I(2)-1;                  
             else 
                 new_row = I(1)+1;
                 new_column = I(2)+1;
             end
         end
    end
end

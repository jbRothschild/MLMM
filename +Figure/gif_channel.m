function gif_channel(t, x, p, colours)
filename = 'trial.gif';
import Figure.*

delaytime = 0.5;
time_gif = 30.;

if size(t,1) >= time_gif
    delaytime = delaytime*time_gif/size(t,1);
    
for i = 1:size(t,1)
    % Capture the plot as an image 
    
    frame = getframe(channel_heatmap(i, x, p, colours, false)); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delaytime); 
    else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delaytime); 
    end 
    close
end

end
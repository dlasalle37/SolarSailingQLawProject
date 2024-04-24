clear; close all

%% Process Data
angles = readmatrix("data\angles.txt"); % sail angles [alpha, beta]
t = readmatrix("data\discrete_times.txt"); % time vector, starting at 0

interval = 7; %minutes
N = length(angles);
set(gca, 'Xlim', [0, t(end)], 'Ylim', [0, 90])
ylabel('Angle [deg]'); xlabel('Time [s]')
curve = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5);
f = cell(N, 1);
k=1;
for i=1:interval:N
    addpoints(curve, t(i), angles(i,1)*180/pi)
    drawnow
    f{k} = getframe(gcf) ;
    k=k+1;
end

obj = VideoWriter('curve.avi');
obj.Quality = 100;
obj.FrameRate = 2;
open(obj);
for i = 1:k-1
    writeVideo(obj, f{i}) ;
end
obj.close();
hold off
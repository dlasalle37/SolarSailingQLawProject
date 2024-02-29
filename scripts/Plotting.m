% Plotting
% - Copy data from Julia Simulations (in highest-level data folder)
%   into the data folder in this directory to plot. This is done so the
%   data here is never accidentally overwritten when trying to compare
%   results vs the Julia Plotter.

clear; close all
r_body = 6378.0; % changed based on central body desired
ecl = readmatrix("data\eclipsed.txt"); % coasting arcs
cart = readmatrix("data\cart.txt"); % non-coasting arcs
kep = readmatrix("data\kep.txt"); % Keplerian elements
angles = readmatrix("data\angles.txt"); % sail angles [alpha, beta]
t = readmatrix("data\discrete_times.txt"); % time vector, starting at 0

sP = [cart(1,1); cart(1,2); cart(1, 3)]; % start point
eP = [cart(end,1); cart(end,2); cart(end, 3)]; % end point
% Figure 1 ==========
figure
hold on
plot3(cart(:,1), cart(:,2), cart(:,3))
plot3(ecl(:,1), ecl(:,2), ecl(:,3), '.r', 'markersize', 3)
plot3(sP(1), sP(2), sP(3), 'xk', 'markersize', 6, 'linewidth', 1)
plot3(eP(1), eP(2), eP(3), 'sk', 'markersize', 6, 'linewidth', 1) 
[x, y, z] = sphere(25); % Making the earth Sphere
surf(r_body*x, r_body*y, r_body*z)
title('Trajectory Plot, Inertial')
xlabel("Inertial X[km]"); ylabel("Inertial Y[km]"); zlabel("Inertial Z[km")
view(30, 30)
axis equal
legend("Trajectory", "Eclipsed Points", "Starting Point", "Ending Point", " ")
hold off

% Figure 2 ==========
figure

subplot(3, 2, 1)
plot(t/86400, kep(:,1))
title("Semi-major axis [km]")
xlabel("t [days]")

subplot(3, 2, 2)
plot(t/86400, kep(:,2))
xlabel("t [days]")
title("Eccentricity")

subplot(3, 2, 3)
plot(t/86400, kep(:,3)*180/pi)
xlabel("t [days]")
title("Inclination [deg]")

subplot(3, 2, 4)
plot(t/86400, kep(:,4)*180/pi)
xlabel("t [days]")
title("Arg. perigee [deg]")

subplot(3, 2, 5)
plot(t/86400, kep(:,5)*180/pi)
xlabel("t [days]")
title("RAAN [deg]")

subplot(3, 2, 6)
plot(t/86400, kep(:,6)*180/pi)
xlabel("t [days]")
title("True anom [deg]")

figure(3)
subplot(2,1,1)
plot(t/86400, angles(:,1)*180/pi)
title("Sail Angles")
xlabel('Time [days]'); ylabel("Cone angle \alpha [degrees]")

subplot(2,1,2)
plot(t/86400, angles(:,2)*180/pi)
xlabel('Time [days]'); ylabel("Clock angle \beta [degrees]")
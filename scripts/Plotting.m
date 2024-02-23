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
t = readmatrix("data\discrete_times.txt"); % time vector, starting at 0

% Figure 1 ==========
figure
hold on
plot3(cart(:,1), cart(:,2), cart(:,3))
plot3(ecl(:,1), ecl(:,2), ecl(:,3), '.r', 'markersize', 3)
[x, y, z] = sphere(25); % Making the earth Sphere
surf(r_body*x, r_body*y, r_body*z)
title('Trajectory Plot, Inertial')
xlabel("Inertial X[km]"); ylabel("Inertial Y[km]"); zlabel("Inertial Z[km")
view(30, 30)
axis equal
legend("Trajectory", "Eclipsed Points", " ")
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

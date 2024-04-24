%%  Plotting
% - Copy data from Julia Simulations (in highest-level data folder)
%   into the data folder in this directory to plot. This is done so the
%   data here is never accidentally overwritten when trying to compare
%   results vs the Julia Plotter.
clear; close all

%% Settings
transferName    = "A";                  % which transfer to plot
font            = "Times New Roman";    %desired font
fontSize        = 11;                   % font size
saveFigs        = true;                % do you want to save output (.eps) figures
height          = 3.25;                 % for saved figures
%% Process Data
r_body = 6378.0; % changed based on central body desired

ecl = readmatrix("data\plotting_data\Case"+transferName+"_eclipsed.txt"); % coasting arcs
cart = readmatrix("data\plotting_data\Case"+transferName+"_cart.txt"); % non-coasting arcs
kep = readmatrix("data\plotting_data\Case"+transferName+"_kep.txt"); % Keplerian elements
angles = readmatrix("data\plotting_data\Case"+transferName+"_angles.txt"); % sail angles [alpha, beta]
t = readmatrix("data\plotting_data\Case"+transferName+"_discrete_times.txt"); % time vector, starting at 0


sP = [cart(1,1); cart(1,2); cart(1, 3)]; % start point
eP = [cart(end,1); cart(end,2); cart(end, 3)]; % end point

% Creating Initial/final orbits
mu = 398600.4418; 
X0 = cart(1,:); XF = cart(end,:);
func = @(t, s) twobeom(t, s, mu);
ai = kep(1,1); af = kep(end, 1);
[T, Xi] = ode45(func, [0 2*pi*sqrt(ai^3/mu)], X0, odeset('reltol', 1.0e-6, 'abstol', 1.0e-6));
[T, Xf] = ode45(func, [0 2*pi*sqrt(af^3/mu)], XF, odeset('reltol', 1.0e-6, 'abstol', 1.0e-6));

%% Make figures

%% Figure 1 ==========
figure(1)
subplot(3, 2, 1)
plot(t/86400, kep(:,1))
ylabel("a [km]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")


subplot(3, 2, 2)
plot(t/86400, kep(:,2))
ylabel("e", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")


subplot(3, 2, 3)
plot(t/86400, kep(:,3)*180/pi)
ylabel("i [deg]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")


subplot(3, 2, 4)
plot(t/86400, kep(:,4)*180/pi)
ylabel("$\omega$ [deg]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")


subplot(3, 2, 5)
plot(t/86400, kep(:,5)*180/pi)
ylabel("$\Omega$ [deg]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
xlabel("t [days]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")


subplot(3, 2, 6)
plot(t/86400, kep(:,6)*180/pi)
xlabel("t [days]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
ylabel("$\nu$ [deg]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);

% Setup for publication figures
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")

if saveFigs==true
    print("Case"+transferName+"_oe",'-depsc2');
end
%%  Figure 2
figure(2)
subplot(2,1,1)
plot(t/86400, angles(:,1)*180/pi)
xlabel("Time [days]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
ylabel("$\alpha$ [degrees]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")

subplot(2,1,2)
plot(t/86400, angles(:,2)*180/pi)
xlabel("Time [days]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
ylabel("$\beta$ [degrees]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")
if saveFigs==true
    print("Case"+transferName+"_angles",'-depsc2');
end
%% Figure 3 =========
figure(3)
hold on
plot3(cart(:,1), cart(:,2), cart(:,3))
plot3(ecl(:,1), ecl(:,2), ecl(:,3), '.r', 'markersize', 3)
% plot3(sP(1), sP(2), sP(3), 'xk', 'markersize', 6, 'linewidth', 1)
% plot3(eP(1), eP(2), eP(3), 'sk', 'markersize', 6, 'linewidth', 1) 
plot3(Xi(:,1), Xi(:,2), Xi(:,3), '-k', 'linewidth', 2)
plot3(Xf(:,1), Xf(:,2), Xf(:,3), '-g', 'linewidth', 2)
[x, y, z] = sphere(25); % Making the earth Sphere
surf(r_body*x, r_body*y, r_body*z, 'HandleVisibility', 'off')
xlabel("Inertial X[km]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font); 
ylabel("Inertial Y[km]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font); 
zlabel("Inertial Z[km]", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
axis equal
legend("Trajectory", "Eclipsed Points",'Initial Orbit',...
    'Final Orbit', " ", "Location", "eastoutside")
% Setup for publication figures
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
set(gcf, "PaperPositionMode","Manual")
hold off
view(0, 90)
if saveFigs==true
    print("Case"+transferName+"_cart",'-depsc2');
end
function dxdt = twobeom(t, s, mu)
    rvec = s(1:3);
    vvec = s(4:6);
    r = norm(rvec);
    
    dxdt = [vvec; (-mu/r^3) * rvec];
end
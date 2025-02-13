%%  Plotting
% - Copy data from Julia Simulations (in highest-level data folder)
%   into the data folder in this directory to plot. This is done so the
%   data here is never accidentally overwritten when trying to compare
%   results vs the Julia Plotter.
clear; close all

%% Settings
transferName    = "C";                  % which transfer to plot
font            = "Times New Roman";    %desired font
fontSize        = 10;                   % font size
saveFigs        = false;                 % do you want to save output (.eps) figures
height          = 2.667;                 % for saved figures
width           = 4;
mdl             = "5x5";                % grav model used (5x5 or 2x0)
j2              = "WITH";                 % "WITH" or "NO"

%% Process Data
r_body = 6378.0; % changed based on central body desired

ecl = readmatrix("data\plotting_data\"+mdl+"\CASE_"+transferName+"_5by5_ECLIPSED_"+j2+"_J2.txt"); % coasting arcs
cart = readmatrix("data\plotting_data\"+mdl+"\CASE_"+transferName+"_5by5_"+j2+"_J2.txt"); % non-coasting arcs
kep = readmatrix("data\plotting_data\Case"+transferName+"_kep.txt"); % Keplerian elements
angles = readmatrix("data\plotting_data\"+mdl+"\CASE_"+transferName+"_5by5_ANGLES_"+j2+"_J2.txt"); % sail angles [alpha, beta]
t = 0:60:(60*length(cart)-1);

kep_no_J2_in_qlaw = readmatrix("data\plotting_data\"+mdl+"\CASE_"+transferName+"_5by5_NO_J2_oe.txt");
kep_with_J2_in_qlaw = readmatrix("data\plotting_data\"+mdl+"\CASE_"+transferName+"_5by5_WITH_J2_oe.txt");
%kep_with_J2_in_qlaw(:,5) = mod(kep_with_J2_in_qlaw(:,5), 2*pi);
t_no_j2 = 0:60:(60*length(kep_no_J2_in_qlaw)-1);
t_with_j2 = 0:60:(60*length(kep_with_J2_in_qlaw)-1);
% Q_with_j2 = readmatrix("data\plotting_data\CASE_"+transferName+"_Q_WITH_J2.txt");
% Q_no_j2 = readmatrix("data\plotting_data\CASE_"+transferName+"_Q_NO_J2.txt");

sP = [cart(1,1); cart(1,2); cart(1, 3)]; % start point
eP = [cart(end,1); cart(end,2); cart(end, 3)]; % end point

% Creating Initial/final orbits
mu = 398600.4418; 
X0 = cart(1,:); XF = cart(end,:);
func = @(t, s) twobeom(t, s, mu);
ai = kep(1,1); af = kep(end, 1);
[T, Xi] = ode45(func, [0 2*pi*sqrt(ai^3/mu)], X0, odeset('reltol', 1.0e-6, 'abstol', 1.0e-6));
[T, Xf] = ode45(func, [0 2*pi*sqrt(af^3/mu)], XF, odeset('reltol', 1.0e-6, 'abstol', 1.0e-6));

%% Pull case info
switch(transferName)
    case("A")
        Woe = [1, 1, 1, 0, 0];
        oet = [26500.0, 0.75, 0.01*pi/180, 0.0*pi/180, 90*pi/180];
        ylim_a = [kep_no_J2_in_qlaw(1,1)-1500, oet(1)+1500];
        ylim_e = [kep_no_J2_in_qlaw(1,2)-0.075, oet(2)+0.075];
        ylim_i = [min(kep_no_J2_in_qlaw(:,3))*180/pi-.1, max(kep_no_J2_in_qlaw(:,3))*180/pi];
    case("B")
        Woe = [1.0, 1.0, 1.0, 1.0, 0.0];
        oet = [26700, 0.75, 0.2*pi/180, 30*pi/180, 90.0*pi/180];
        ylim_a = [kep_no_J2_in_qlaw(1,1)-1500, oet(1)+1500];
        ylim_e = [kep_no_J2_in_qlaw(1,2)-0.075, oet(2)+0.075];
        ylim_i = false;
    case("C")
        oet = [26700, 0.205, 60.0*pi/180, 30*pi/180, 30*pi/180];
        Woe = [10.0, 10.0, 1.0, 1.0, 5.0];
        ylim_a = [25000, 30000];
        ylim_e = [0, 0.5];
        ylim_i = false;
end

%% Figure 1 ==========
figure(1)
tcl = tiledlayout(3,2);
nexttile(tcl)
hold on; box on
line1 = plot(t_no_j2/86400, kep_no_J2_in_qlaw(:,1), 'DisplayName', 'Without J2');
line2 = plot(t_with_j2/86400, kep_with_J2_in_qlaw(:,1), 'DisplayName', 'With J2');
pvec = ones(length(kep_no_J2_in_qlaw), 1);
if Woe(1) ~= 0
    line3 = plot(t_no_j2/86400, pvec*oet(1), '--k', 'DisplayName', 'Target');
end
ylabel("a, km", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25, 5, 3.75])
set(gcf, "PaperPositionMode","Manual")
ylim(ylim_a)
hold off

nexttile(tcl)
hold on; box on
plot(t_no_j2/86400, kep_no_J2_in_qlaw(:,2), t_with_j2/86400, kep_with_J2_in_qlaw(:,2))
if Woe(2) ~= 0
    plot(t_no_j2/86400, pvec*oet(2), '--k')
end
ylabel("e", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,5,3.75])
set(gcf, "PaperPositionMode","Manual")
ylim(ylim_e)
hold off

nexttile(tcl)
hold on; box on
plot(t_no_j2/86400, kep_no_J2_in_qlaw(:,3)*180/pi, t_with_j2/86400, kep_with_J2_in_qlaw(:,3)*180/pi)
if Woe(3) ~= 0
    plot(t_no_j2/86400, pvec*oet(3)*180/pi, '--k')
end
ylabel("i, degrees", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
if ylim_i ~= false
    ylim(ylim_i)
end
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25, 5, 3.75])
set(gcf, "PaperPositionMode","Manual")
hold off

nexttile(tcl)
hold on; box on
plot(t_no_j2/86400, kep_no_J2_in_qlaw(:,4)*180/pi, t_with_j2/86400, kep_with_J2_in_qlaw(:,4)*180/pi)
if Woe(4) ~= 0
    plot(t_no_j2/86400, pvec*oet(4)*180/pi, '--k')
end
ylabel("$\omega$, degrees", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25, 5, 3.75])
set(gcf, "PaperPositionMode","Manual")
%ylim([-60 35])
hold off

nexttile(tcl)
hold on; box on
plot(t_no_j2/86400, kep_no_J2_in_qlaw(:,5)*180/pi, t_with_j2/86400, kep_with_J2_in_qlaw(:,5)*180/pi)
if Woe(5) ~= 0
    plot(t_no_j2/86400, pvec*oet(5)*180/pi, '--k')
end
%ylim([-10 35])
ylabel("$\Omega$, degrees", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
xlabel("t, days", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25, 5, 3.75])
set(gcf, "PaperPositionMode","Manual")
hold off

nexttile(tcl)
plot(t_no_j2/86400, kep_no_J2_in_qlaw(:,6)*180/pi, t_with_j2/86400, kep_with_J2_in_qlaw(:,6)*180/pi)
xlabel("t, days", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
ylabel("$\nu$, degrees", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);

lg = legend([line1,line2,line3], 'Location', 'northoutside');
lg.Orientation = 'horizontal';
% Setup for publication figures
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25, 5, 3.75])
set(gcf, "PaperPositionMode","Manual")

if saveFigs==true
    print("Case"+transferName+"_oe",'-depsc2');
end
%%  Figure 2
figure(2)
subplot(2,1,1)
plot(t/86400, angles(:,1)*180/pi)
xlabel("t, days", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
ylabel("$\alpha$, degrees", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,2.5,2])
set(gcf, "PaperPositionMode","Manual")

subplot(2,1,2)
plot(t/86400, angles(:,2)*180/pi)
xlabel("t, days", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
ylabel("$\beta$, degrees", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,2.5,2])
set(gcf, "PaperPositionMode","Manual")
if saveFigs==true
    print("Case"+transferName+"_angles_"+j2+"_j2",'-depsc2');
end
%% Figure 3 =========
figure(3)
hold on

% Earth plotting
[xS,yS,zS] = sphere(50);
r_earth = 6378.137; %km
xSE = r_earth*xS;
ySE = r_earth*yS;
zSE = r_earth*zS;
s = surface(xSE,ySE,zSE); % Plot the Earth Shere
axis equal %axes must be set to equal
%Wrapping the map:
imData = imread('2_no_clouds_4k.jpg');
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',flipud(imData),'edgecolor','none');
s.HandleVisibility = 'off'; % turn off surface handle vis so it doesnt show in legend

plot3(cart(:,1), cart(:,2), cart(:,3), 'linestyle', '-', 'color', '#4169e1')
plot3(ecl(:,1), ecl(:,2), ecl(:,3), '.r', 'markersize', 3)
% plot3(sP(1), sP(2), sP(3), 'xk', 'markersize', 6, 'linewidth', 1)
% plot3(eP(1), eP(2), eP(3), 'sk', 'markersize', 6, 'linewidth', 1) 
plot3(Xi(:,1), Xi(:,2), Xi(:,3), '-k', 'linewidth', 2)
plot3(Xf(:,1), Xf(:,2), Xf(:,3), '-g', 'linewidth', 2)
% [x, y, z] = sphere(25); % Making the earth Sphere
% surf(r_body*x, r_body*y, r_body*z, 'HandleVisibility', 'off')
xlabel("Inertial X, km", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font); 
ylabel("Inertial Y, km", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font); 
zlabel("Inertial Z, km", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
axis equal; box on

% Custom legend
cm{1} = plot(nan, 'linestyle', '-', 'color', '#4169e1', 'linewidth', 1.5);
cm{2} = plot(nan, '-r', 'linewidth', 1.5);
cm{3} = plot(nan, '-k', 'linewidth', 1.5);
cm{4} = plot(nan, '-g', 'linewidth', 1.5);
legend([cm{:}], {'Thrusting', 'Eclipsed', 'Initial Orbit', 'Final Orbit'},...
    'Location', 'eastoutside', 'fontsize', fontSize)
% legend("Trajectory", "Eclipsed Points",'Initial Orbit',...
%         'Final Orbit', " ", "Location", "eastoutside", "fontsize", fontSize)
    
% Setup for publication figures
set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25, 5, height])
% set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,2.5,2])
set(gcf, "PaperPositionMode","Manual")
view(0, 90)
view(90, 0)

if saveFigs==true
    print("Case"+transferName+"_cart_"+j2+"_j2",'-depsc2', '-r300');
end

% %% Figure 4
% figure(4)
% plot(t_no_j2/86400, Q_no_j2, t_with_j2/86400, Q_with_j2)
% xlabel("Time (days)")
% ylabel("Q")
% % Setup for publication figures
% set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
% set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.68504,height])
% set(gcf, "PaperPositionMode","Manual")
% % if saveFigs==true
% %     print("Case"+transferName+"_Q",'-depsc2');
% % end 

% %% figure 4
% figure
% hold on; box on
% line1 = plot(t_no_j2/86400, kep_no_J2_in_qlaw(:,1), 'DisplayName', 'Without J2');
% line2 = plot(t_with_j2/86400, kep_with_J2_in_qlaw(:,1), 'DisplayName', 'With J2');
% pvec = ones(length(kep_no_J2_in_qlaw), 1);
% if Woe(1) ~= 0
%     line3 = plot(t_no_j2/86400, pvec*oet(1), '--k', 'DisplayName', 'Target');
% end
% ylabel("a, km", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
% set(gca, "fontname", "Times New Roman", "fontsize", fontSize)
% set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,4.5,height])
% set(gcf, "PaperPositionMode","Manual")
% ylim(ylim_a)
% hold off
% xlim([20, 70])
% ylim([2.4e4, 2.8e4])
% xlabel("t, days", "Interpreter", "Latex", "FontSize", fontSize, "fontname", font);
% legend("Without J2", "With J2", "Target")
% %view(90, 0)
% if saveFigs==true
%     print("Case"+transferName+"_oe_zoomed", '-depsc2');
% end

function dxdt = twobeom(t, s, mu)
    rvec = s(1:3);
    vvec = s(4:6);
    r = norm(rvec);
    
    dxdt = [vvec; (-mu/r^3) * rvec];
end
clear; close all

%% Process Data
r_body = 6378.0; % changed based on central body desired
ecl = readmatrix("data\eclipsed.txt"); % coasting arcs
cart = readmatrix("data\cart.txt"); % non-coasting arcs
kep = readmatrix("data\kep.txt"); % Keplerian elements
angles = readmatrix("data\angles.txt"); % sail angles [alpha, beta]
t = readmatrix("data\discrete_times.txt"); % time vector, starting at 0

sP = [cart(1,1); cart(1,2); cart(1, 3)]; % start point
eP = [cart(end,1); cart(end,2); cart(end, 3)]; % end point

% Creating Initial/final orbits
mu = 398600.4418; 
X0 = cart(1,:); XF = cart(end,:);
func = @(t, s) twobeom(t, s, mu);
ai = kep(1,1); af = kep(end, 1);
[T, Xi] = ode45(func, [0 2*pi*sqrt(ai^3/mu)], X0, odeset('reltol', 1.0e-6, 'abstol', 1.0e-6));
[T, Xf] = ode45(func, [0 2*pi*sqrt(af^3/mu)], XF, odeset('reltol', 1.0e-6, 'abstol', 1.0e-6));

%% Make animation
figure
hold on
X = cart(:,1); Y = cart(:,2); Z = cart(:,3);
N = length(X);
interval = 7; % minutes
curve = animatedline('Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
f = cell(length(X), 1);
set(gca)
xlabel('Inertial X [km]'); ylabel('Inertial Y [km]');
plot3(Xf(:,1), Xf(:,2), Xf(:,3), '-g', 'linewidth', 2)
plot3(Xi(:,1), Xi(:,2), Xi(:,3), '-k', 'linewidth', 2)
[x, y, z] = sphere(25); % Making the earth Sphere
surf(r_body*x, r_body*y, r_body*z)
axis equal
k=1;
for i=1:interval:N
    addpoints(curve, X(i), Y(i), Z(i))
    head = scatter3(X(i), Y(i), Z(i), "+r", "linewidth", 2);
    drawnow
    f{k} = getframe(gcf) ;
    delete(head);
    k=k+1;
end
legend("Initial", "Final")
obj = VideoWriter('curve.avi');
obj.Quality = 100;
obj.FrameRate = 2;
open(obj);
for i = 1:k-1
    writeVideo(obj, f{i}) ;
end
obj.close();
hold off
%% Supplementary
function dxdt = twobeom(t, s, mu)
    rvec = s(1:3);
    vvec = s(4:6);
    r = norm(rvec);
    
    dxdt = [vvec; (-mu/r^3) * rvec];
end
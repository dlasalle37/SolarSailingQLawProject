% two-body orbit prop
clear; close all
mu = 398600.4418;
coe0 = [10500.0, 0.150, 25.0, 30.0, 10.0, 0.0];
[r, v] = coe2rv(coe0, mu);
X0 = [r;v];
tspan = 0:60:86400*5; 
odefun = @(t, s) twobeom(t, s, mu);
[T, X] = ode45(odefun, tspan, X0, odeset('reltol', 1.0e-12, 'abstol', 1.0e-12)); 

writematrix(X,'unit_test_scripts\data_test\MatlabSimData.txt','Delimiter','\t')  

%% Functions
function dxdt = twobeom(t, s, mu)
    rvec = s(1:3);
    vvec = s(4:6);
    r = norm(rvec);
    
    dxdt = [vvec; (-mu/r^3) * rvec];
end

function coe_set = rv2coe(r, v,mu)

    r_norm = norm(r); %scalar position
    v_norm = norm(v); %scalar velocity

    n_z = [0 0 1]'; %inertial z unit vector
    n_x = [1 0 0]'; %inertial x unit vector

    %Semi-Major
    %From vis-viva:
    a = (2/r_norm - v_norm^2/mu)^-1;

    %Eccentricity
    h_vec = cross(r,v); %angular momentum vector
    e_vec = 1/mu * (cross(v,h_vec)) - 1/r_norm * r; %eccentricity vector
    e = norm(e_vec); %eccentricity

    %true anomaly
    nu = acosd(dot(r, e_vec) / (r_norm * e)); % true anomaly
    %need to check condition
    %if dot(r,v)>0, 0<nu<180 (no change needed)
    %if dot(r,v)<0, 180<nu<360

    if dot(r,v)<0
        nu = 360-nu;
    end

    %Inclination
    h3  = h_vec(3); %z component of h
    h = norm(h_vec); %scalar ang. momentum
    i = acosd(h3/h); 

    %RAAN
    h = norm(h_vec);
    e_h = h_vec / h; %e_h vector from eq (115) 
    e_n = cross(n_z, e_h);
    en1 = dot(e_n, n_x); %grab first (x) component of e_n
    RAAN = acosd(en1/norm(e_n)); 

    %conditional for second component of e_n being negative
    if e_n(2) < 0
        RAAN = 360 - RAAN;
    end

    %argument of perigee
    e_e = e_vec / e; %normalized eccentricity vector in inertial coords
    w = acosd( dot(e_n,e_e) / (norm(e_n) * norm(e_e)) );
    %conditional:
    if e_e(3) < 0
        w = 360 - w;
    end

    coe_set = [a, e, i, RAAN, w, nu];
end

function [r_n, v_n] = coe2rv(coe_set, mu)
% Unpack COE Set
a = coe_set(1); e = coe_set(2); i = coe_set(3);
RAAN = coe_set(4); w = coe_set(5); nu = coe_set(6);

p = a*(1-e^2); %semi-latus rectum
r = p / (1+e*cosd(nu)); %scalar r from trajectory

r_p = [r*cosd(nu); r*sind(nu); 0]; %perifocal position vector


v_p = [-sqrt(mu/p) * sind(nu); sqrt(mu/p)*(e+cosd(nu)); 0]; %perifocal v

%need to transform r_p and v_p  to inertial
%use p->n transfer matrix, A. (eq126):

A11 = cosd(RAAN)*cosd(w) - sind(RAAN)*sind(w)*cosd(i);
A12 = -cosd(RAAN)*sind(w) - sind(RAAN)* cosd(w)*cosd(i);
A13 = sind(RAAN)*sind(i);
A21 = sind(RAAN)*cosd(w) + cosd(RAAN)*sind(w)*cosd(i);
A22 = -sind(RAAN)*sind(w) + cosd(RAAN)*cosd(w)*cosd(i);
A23 = -cosd(RAAN)*sind(i);
A31 = sind(w)*sind(i);
A32 = cosd(w)*sind(i);
A33 = cosd(i);

A = [A11 A12 A13;
     A21 A22 A23;
     A31 A32 A33];
%inertial vectors:

r_n = A*r_p; 
v_n = A*v_p;
 
end




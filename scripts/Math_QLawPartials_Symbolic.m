%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the derivatives of Q with respect to the 5 slow elements:
% semi-major axis, eccentricity,  inclination, argument of perigee, and 
% lambda = RAAN-trueAnom_Earth
% Symbolic calculation, formulated from "Solar Sailing Q-Law for 
% Planetocentric, Many-Revolution Sail Orbit Transfers" (Oguri et. al)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definitions
clear
% Define Symbolic variables
syms Q dQdt real  % Difference quotient and time derivative
syms mu mu_sun real  % gravitational parameter of central bodies
syms e_E a_E tru_E real% earth's helio coe's
syms a_SRP real         % magnitude of control acceleration
syms a_SRP_r a_SRP_theta a_SRP_h real  % components of control acceleration in orbit frame
syms alpha beta real  % control angles alpha and beta
syms h p r real % osculating values of: angular momentum magnitude, semi-latus rectum, and position vector magnitude
syms syms W Wa We Winc Wape Wlam real  %COE selection weights

syms P Aimp Aesc kimp kesc Wp rp rpmin amax real % Penalty function stuff

syms a e inc ape lam tru real         % osculating orbital elements (and true anom)

syms a_t e_t inc_t ape_t lam_t real   % target orbital elements

syms Qa Qe Qinc Qape Qlam real        % terms of Lyapunov function summation
syms m_petro n_petro real
syms r_petro b_petro k_petro real
syms C1 C2 C3 Asail msail real % sail parameters
syms r_sc_sun G0 real% SRP parameters (assuming r_sc_sun is approximately constant)
syms nustar_a nustar_e nustar_inc nustar_ape nustar_lam real
% ^ terms to calculate oe_dot_xx

syms dQdcoe real% 5x5 array holding derivatives of each term of

% Element selection vectors
eps_a = [1 0 0 0 0 0]';
eps_e = [0 1 0 0 0 0]'; 
eps_inc = [0 0 1 0 0 0]';
eps_ape = [0 0 0 1 0 0]';
eps_lam = [0 0 0 0 1 0]';

%% Lyapunov function basics
% Lyapunov function with respect to each orbital element
% put h, p, r in terms of orbital elements
p = a*(1-e^2);
h = sqrt(mu*p);
r = p/(1+e*cos(tru));
rp = a*(1-e);

% orbit to hill transformation matrix
R_H_O = rot_hill_to_orbit(inc, ape, lam, tru);


% Scaling function (semi-major axis only)
Sa = (1+((a-a_t)/(m_petro*a_t))^n_petro)^(1/r_petro);

% Penalty function
cimp = 1 - a*(1-e)/rpmin;  % impact constraint
cesc = a/amax - 1;  % escape constraint
P = Aimp*exp(kimp*cimp) + Aesc*exp(kesc*cesc);

%%  Calculating parameters for best-case time-to-go
% Distance Functions
dista = a - a_t;
diste = e - e_t;
distinc = inc - inc_t;
distape = acos(cos(ape - ape_t));
distlam = acos(cos(lam-lam_t));

% sigma functions 
sig_a = sign(dista);
sig_e = sign(diste);
sig_inc = sign(distinc);
sig_ape = sign(distape);
sig_lam = sign(distlam);

% Sunlight direction in perifocal frame
term1 = [cos(ape) sin(ape) 0; -sin(ape) cos(ape) 0; 0 0 1];
term2 = [1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)];
term3 = [cos(lam) sin(lam) 0; -sin(lam) cos(lam) 0; 0 0 1];
svec_P = term1*term2*term3*[1;0;0];
se = svec_P(1); sp = svec_P(2); sh = svec_P(3); %breaking it down

% ballistic evolution of state
nudot = (1+e*cos(tru))^2 / (1-e^2)^(3/2) * sqrt(mu/a^3);
nudot_body = (1+e_E*cos(tru_E))^2 / (1-e_E^2)^(3/2) * sqrt(mu_sun/a_E^3);
f0 = [0 0 0 0 -nudot_body nudot]'; %f0(xslow)

% Calculating oedot_nn for each oe
% Semi-major axis:
R_H_O_star_a = rot_hill_to_orbit(inc, ape, lam, nustar_a);
F_a = F(a, e, inc, ape, lam, nustar_a, mu);
p_a = -[sig_a*eps_a'*F_a*R_H_O_star_a]';
p_ax = p_a(1); p_ay = p_a(2); p_az = p_a(3);

betastar_a = atan2(-p_ay, -p_az);
k = p_ax/sqrt(p_ay^2 + p_az^2);
alphstar_a0 = atan(1/4 * (-3*k + sqrt(8+9*k^2)));
F_alpha_a = k*(3*C1*cos(alphstar_a0)^2 + 2*C2*cos(alphstar_a0) + C3)*sin(alphstar_a0)...
    -C1*cos(alphstar_a0)*(1-3*sin(alphstar_a0)^2) - C2*cos(2*alphstar_a0);
F_alpha_alpha_a = k*(3*C1*(cos(alphstar_a0)-2*cos(alphstar_a0))*sin(alphstar_a0)...
    + 2*C2*cos(2*alphstar_a0) + C3*cos(alphstar_a0))...
    -C1*sin(alphstar_a0)*(2-9*cos(alphstar_a0))+2*C2*sin(2*alphstar_a0);
alphastar_a = alphstar_a0 - F_alpha_a/F_alpha_alpha_a;
ustar_a = [alphastar_a; betastar_a];
astar_SRP_a = a_SRP_hill(ustar_a, Asail, msail, G0, r_sc_sun, C1, C2, C3);
adotnn = sig_a*eps_a'*f0-p_a'*astar_SRP_a;
Qa = (dista^2 / (-abs(dista)*adotnn))^2;

% Eccentricity: use both, see which gives greater 

% Inclination

% Lam


% Argument of periapsis:
% So far I have found two ways of doing this as there is not yet an
% analytical solution:
% 1.) solve for in-plane and out of plane motions of ape_dot, weight with b
% 2.) numerically solve for the nu that minimizes the below equation

% nustar_ape = argmin(...

%% Create Q and dQ/dx

%Delete these as each Qoe is completed
Qe = 0;
Qinc = 0;
Qape=0;
Qlam=0;

Q = (1 + Wp*P)*(Wa*Qa + We*Qe + Winc*Qinc + Wape*Qape * Wlam*Qlam);
dQda = diff(Q,a);

%% Functions
function Fout = F(a, e, inc, ape, lam, tru, mu)
    % F function defined by gauss's variational equations
    
    % Calculate parameters:
    p = a*(1-e^2);
    h = sqrt(mu*p);
    r = p/(1+e*cos(tru));
    
    % Populate matrix
    Fout = (1/h)* [2*a^2*e*sin(tru) 2*a^2*p/r 0;
    p*sin(tru) (p+r)*cos(tru)+r*e 0;
    0 0 r*cos(tru+ape);
    -p*cos(tru)/e (p+r)*sin(tru)/e -r*sin(tru+ape)/tan(inc);
    0 0 r*sin(tru+ape)/sin(inc);
    p*cos(tru)/e -(p+r)*sin(tru)/e 0];
end

function R = rot_hill_to_orbit(inc, ape, lam, tru)
    % Output: Orbit frame to hill frame transformation matrix
    % inputs:
    % inc: inclination
    % ape: argument of perigee
    % lam: orbital element lambda
    % tru: true anomaly
    term1 = [cos(tru+ape) sin(tru+ape) 0;
         -sin(tru+ape) cos(tru+ape) 0;
         0 0 1];
    term2 = [1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)];
    term3 = [cos(lam) sin(lam) 0; -sin(lam) cos(lam) 0; 0 0 1];
    R = term1*term2*term3;
end

function aSRP = a_SRP_hill(u, A, m, G0, d, C1, C2, C3)
    a = u(1); %alpha
    b = u(2); %beta
    
    aSRP = (A/m)*G0/d^2*cos(a)*[C1*cos(a)^2 + C2*cos(a) + C3;
        -(C1*cos(a)+C2)*sin(a)*sin(b);
        -(C1*cos(a)+C2)*sin(a)*cos(b)];
end
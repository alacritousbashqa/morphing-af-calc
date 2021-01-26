% Script file for running Thwaites_Panel function on a geometry
clear variables
close all
clc

N = 200;                % Number of panels
foil = '0012';          % Airfoil four digit NACA number

N = N - mod(N,2);       % Ensures that the number of panels is even

% Generate airfoil
[X,Y] = generateNACA4(foil,N);
X = flip(X);
Y = flip(Y);
Y(1) = 0;
Y(end) = 0;

% Scale coordinates by chord
c = 1;  % Chord [m]
X = c*X;
Y = c*Y;
b = 1;  % Span [m]

rho = 1.225;                    % Freestream Density [kg/m^3]
mu = 1.802e-5;                  % Dynamic Viscosity [kg/m/s]
nu = mu/rho;                    % Kinematic Viscosity [m^2/s]
Vinf = 100;                     % Freestream velocity [m/s]
Re = rho*Vinf*c/mu;             % Chord Reynolds Number

% Angles of attack to run at
alphas = 0;
Cls = zeros(length(alphas),1);
Cds = Cls;

% Run for multiple angles of attack
for k=1:length(alphas)

alpha = alphas(k);
fprintf("Alpha:      %0.2f deg\n", alpha)

% Run a viscous panel method at this AoA
[d_star_t_thw,d_star_b_thw,I_crit,Cp,Cf_t,strengths,X_mid,Y_mid,norms,extra] = Thwaites_panel_1(X,Y,Vinf,alpha,mu,rho,1,false,false);

% Skin friction coefficient
cf_b = Cf_t{1};
cf_t = Cf_t{2};

Cf = [flip(cf_b); cf_t];
Cf(1) = 0;
Cf(isinf(Cf)) = 0;


%% Calculate forces
% Number of points to use on the trailing edge...
top_n = N/2; %... on the top
bot_n = N/2 + 1; %... on the bottom

F_dist = zeros(top_n+bot_n-2,2);    % Force distribution [N]
x_mids = zeros(top_n+bot_n-2,1);    % Panel centers in region of interest (x)
y_mids = x_mids;                    % Panel centers "   "   "   "   "   " (y)
CPs = x_mids;                       % Pressure coefficients "   "   "   "   "   "
dxs = x_mids;                       % Delta x "   "   "   "   "   "
dys = x_mids;                       % Delta y "   "   "   "   "   "
iii = 1;                            % Iteration variable

% x and y lengths of every panel
dx = abs(X(2:end) - X(1:end-1));
dy = abs(Y(2:end) - Y(1:end-1));

% Get the slope of every panel
dy_dx = zeros(N,1);
for j = 1:N
    dy_dx(j) = (Y(j+1)-Y(j))/(X(j+1)-X(j));
end

% Calculate aerodynamic coefficients at this AoA
Ca = 0;         % Normal coefficient
Cn = 0;         % Axial coefficient
for j = 1:length(Cp)/2
    Ca = Ca + (Cp(N-j+1)*dy_dx(N-j+1) - Cp(j)*dy_dx(j))*dx(j) + (Cf(j) + Cf(N-j+1))*dx(j);
    Cn = Cn + (Cp(j) - Cp(N-j+1))*dx(j) + (Cf(N-j+1)*dy_dx(N-j+1) + Cf(j)*dy_dx(j))*dx(j);
end
Cn = Cn/c;
Ca = Ca/c;

Cm = Cp'*(dx'.*X_mid + dy'.*Y_mid);         % Viscous moment coefficient
Cd = Ca*cosd(alpha) + Cn*sind(alpha);       % Viscous drag coefficient
Cl = Cn*cosd(alpha) - Ca*sind(alpha);       % Viscous lift coefficient

Cls(k) = Cl;
Cds(k) = Cd;

% Dynamic pressure [Pa]
q_inf = 0.5*rho*Vinf^2;

end

% Report the chord Reynolds number and lift coefficient
fprintf("Reynolds Number: %0.0f\n", Re)
fprintf("\nCl: %5.2f\n", Cl)


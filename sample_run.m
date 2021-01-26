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

alpha = 5;      % Angle of attack [deg]
fprintf("Alpha:      %0.2f deg\n", alpha)

% Run a viscous panel method at this AoA
[d_star_t_thw,d_star_b_thw,I_crit,Cp,Cf_t,strengths,X_mid,Y_mid,norms,extra] = Thwaites_panel_1(X,Y,Vinf,alpha,mu,rho,1,false,false);

% Plot pressure coefficient
figure
plot(X_mid./c,Cp)
ylabel('Cp')
xlabel('x/c')
xlim([0 1])
set(gca,'YDir','reverse')

% Report chord Reynolds number
fprintf("Reynolds Number: %0.0f\n", Re)


% Script file for performing the double wake method
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

alpha = 15;      % Angle of attack [deg]
fprintf("Alpha:      %0.2f deg\n", alpha)

% Run a viscous panel method at this AoA
[~,~,I_crit] = Thwaites_panel_1(X,Y,Vinf,alpha,mu,rho,1,false,false);

% Report chord Reynolds number
fprintf("Reynolds Number: %0.0f\n", Re)





function [Cp,strengths,X_mid,Y_mid,norms,tans] = inv_wake(X,Y,Vinf,alpha,extra,I)

M = length(X);          % Number of points
N = M-1;                % Number of panels

Qinf = Vinf*[cosd(alpha) sind(alpha)];  % Freestream velocity vector

[X_mid,Y_mid,norms,tans] = geo_decomp(X,Y); % Acquire geoemtric properties from airfoil

thetas = zeros(N,1);    % Holds all panel orientations with respect to the x-axis

%% Calculate influence coefficients

A = zeros(N+1,N+1);     % Influence coefficient matrix
B = A;
b = zeros(N+1,1);       % Freestream contribution

for i = 1:N
    % Coordinates of ith panel end and control point
    %   Control Points: (x,y)
    %   First End Point: (xi,yi)
    %   Second End Point: (xii,yii)
    xi = X(i);
    xii = X(i+1);
    yi = Y(i);
    yii = Y(i+1);
    x = (xi+xii)/2;                 % Control point x-coordinate
    y = (yi+yii)/2;                 % Control point y-coordinate
    
    ti = atan2(yii-yi,xii-xi);      % Panel orientation angle [rad]
    thetas(i) = ti;     % Store orientations
    
    u2 = 0;
    v2 = 0;
    
    for j = 1:N
        % Coordinates of jth panel end points
        %   First End Point: (xj,yj)
        %   Second End Point: (xjj,yjj)
        xj = X(j);
        xjj = X(j+1);
        yj = Y(j);
        yjj = Y(j+1);
        
        u2_old = u2;
        v2_old = v2;
        
        % Get induced velocity at control point due to jth panel endpoints
        % with vortex strength of 1
        [u1,v1,u2,v2] = getVi(1,1,x,y,xj,yj,xjj,yjj,j==i);
        
        if j==1
            A(i,1) = dot([u1 v1],norms(i,:));
            B(i,1) = dot([u1 v1],tans(i,:));
        elseif j==N
            A(i,N) = dot([u1 v1] + [u2_old v2_old],norms(i,:));
            A(i,M) = dot([u2 v2],norms(i,:));
            B(i,N) = dot([u1 v1] + [u2_old v2_old],tans(i,:));
            B(i,M) = dot([u2 v2],tans(i,:));
        else
            A(i,j) = dot([u1 v1] + [u2_old v2_old],norms(i,:));
            B(i,j) = dot([u1 v1] + [u2_old v2_old],tans(i,:));
        end
    end
    
    % Freestream flow contribution
    b(i) = dot(-Qinf,norms(i,:)) - extra(i);
end

% Kutta Condition + wake conditions
A(M,1) = 1;
A(N,:) = zeros(1,M);
A(N,M) = 1;
A(N,I) = 1;
b(I) = 0;
b(M) = 0;

%% Solve linear system for strengths
strengths = A\b;

%% Get Cp and tangential velocities

Cp = zeros(N,1);
vtan = zeros(N,1);
for i=1:N
    for j=1:M
        vtan(i) = vtan(i) + B(i,j)*strengths(j);
    end
    vtan(i) = vtan(i) + dot(Qinf,tans(i,:));
    Cp(i) = 1-(vtan(i)/Vinf)^2;
end

end
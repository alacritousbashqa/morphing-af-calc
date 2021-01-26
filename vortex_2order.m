% Panel Method - Linearly varying vortices
%
% Thus, we need to solve for N+1 vortex strengths using N+1 equations.
%
% Inputted geometry must start at the trailing edge, wraparound clockwise,
% then end at the trailing edge.

% X:            X coordinates of the geometry nodes
% Y:            Y coordinates of the geometry nodes
% Vinf:         Fluid Freestream velocity
% alpha:        Angle of attack of geometry to freestream
% extra:        Factor due to BL instead of directly displacing the airfoil "boundary"
function [Cp,strengths,X_mid,Y_mid,norms,tans] = vortex_2order(X,Y,Vinf,alpha,extra,I)

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

% Kutta Condition
A(M,1) = 1;
A(M,M) = 1;
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

% AUXILIARY FUNCTIONS

% Get induced velocity from vortices at a point (x,y) due to a panel with 
% end points (xj,yj) and (xjj,yjj), and strengths gj and gjj
function [u1,v1,u2,v2] = getVi(gj,gjj,x,y,xj,yj,xjj,yjj,same)
a = atan2(yjj-yj,xjj-xj);   % Panel orientation angle [rad]

% Transformation Matrix (global to panel)
A = [cos(a) sin(a);
    -sin(a) cos(a)];

% Transform to panel coordinates
tmp = A*[x-xj;y-yj];
xp = tmp(1);
yp = tmp(2);
S = (xjj-xj)*cos(a) + (yjj-yj)*sin(a);

% Distances and angles from enpoints j and jj to point (x,y)
rj = sqrt(xp^2 + yp^2);
rjj = sqrt((xp-S)^2 + yp^2);
tj = atan2(yp,xp);
tjj = atan2(yp,xp-S);

% If the control point is on the same panel as the endpoints (e.g. i==j)
if same==true
    u1 = (-0.5*(xp-S)/S)*gj;
    u2 = (0.5*xp/S)*gjj;
    v1 = (-1/2/pi)*gj;
    v2 = (1/2/pi)*gjj;
else
    u1 = (-(yp*log(rjj/rj) + xp*(tjj-tj) - S*(tjj-tj))/(2*pi*S))*gj;
    u2 = ((yp*log(rjj/rj) + xp*(tjj-tj))/(2*pi*S))*gjj;
    v1 = (-((S-yp*(tjj-tj)) - xp*log(rj/rjj) + S*log(rj/rjj))/(2*pi*S))*gj;
    v2 = (((S-yp*(tjj-tj)) - xp*log(rj/rjj))/(2*pi*S))*gjj;
end

% Convert back to global coordinates
tmp = A'*[u1;v1];
u1 = tmp(1);
v1 = tmp(2);

tmp = A'*[u2;v2];
u2 = tmp(1);
v2 = tmp(2);
end



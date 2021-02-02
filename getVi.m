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
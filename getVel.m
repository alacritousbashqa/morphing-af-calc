%% Gets velocity at a point induced by the freestream and all panels
function [V,u,v] = getVel(x,y,X,Y,strengths,Vinf,alpha)

u = 0;
v = 0;
N = length(X)-1;

for j = 1:N
    % Coordinates of jth panel end points
    % First End Point: (xj,yj)
    % Second End Point: (xjj,yjj)
    xj = X(j);
    xjj = X(j+1);
    yj = Y(j);
    yjj = Y(j+1);
    
    % Get velocity at (x,y) due to panel with nodes (xj,yj) and (xjj,yjj)
    [u1,v1,u2,v2] = getVi(strengths(j),strengths(j+1),x,y,xj,yj,xjj,yjj,(x==(xj+xjj)/2 && y==(yj+yjj)/2));
    
    if j==1
        uu = u1;
        vv = v1;
        uu_old = u2;
        vv_old = v2;
    elseif j==N
        uu = uu_old + u1 + u2;
        vv = vv_old + v1 + v2;
    else
        uu = uu_old + u1;
        vv = vv_old + v1;
        uu_old = u2;
        vv_old = v2;
    end
    
    u = u + uu;
    v = v + vv;
end

% Convert and add freestream component
u = u + Vinf*cosd(alpha);
v = v + Vinf*sind(alpha);
% Calculate velocity magnitude
V = sqrt(u^2 + v^2);
end
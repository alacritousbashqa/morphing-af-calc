%% Gets the control points, normals, and tangents of each panel from the inputted geometry

function [X_mid,Y_mid,norms,tans,node_norms] = geo_decomp(X,Y)
N = length(X)-1;        % Number of panels

X_mid = zeros(N,1);             % Holds all control point X locations
Y_mid = zeros(N,1);             % Holds all control point Y locations
norms = zeros(N,2);             % Holds all panel normals
node_norms = zeros(N+1,2);      % Holds all node normals
tans = norms;                   % Holds all panel tangents

for i = 1:N
    % Coordinates of ith panel end and control points
    %   Control Points: (x,y)
    %   First End Point: (xi,yi)
    %   Second End Point: (xii,yii)
    xi = X(i);
    xii = X(i+1);
    x = (xi+xii)/2;
    X_mid(i) = x;
    yi = Y(i);
    yii = Y(i+1);
    y = (yi+yii)/2;
    Y_mid(i) = y;
    
    ti = atan2(yii-yi,xii-xi);          % Angle between x-axis and panel
    norms(i,:) = [-sin(ti) cos(ti)];    % ith Panel normal
    tans(i,:) = [cos(ti) sin(ti)];      % ith Panel tangent
    % Calculate node normals as averages of surrounding panel normals
    if i==1
        node_norms(1,:) = norms(1,:);
    elseif i==N
        node_norms(i,:) = (norms(i,:)+norms(i-1,:))/2;
        node_norms(N+1,:) = norms(N,:);
    else
        node_norms(i,:) = (norms(i,:)+norms(i-1,:))/2;
    end
end
end
% This file adds in viscous effects by modelling a boundary layer
% across the surfaces of the airfoil based on Thwaites method.
%
% This file using the new linearly varying vortex method for the inviscid 
% solution. This file also includes the geo_decom function for conciseness.

% X:            X coordinates of the geometry nodes
% Y:            Y coordinates of the geometry nodes
% Vinf:         Fluid Freestream velocity
% alpha:        Angle of attack of geometry to freestream
% mu:           Fluid Dynamic viscosity
% rho:          Fluid Density
% it:           Number of iterations to complete
% b_plotFoil:   Show airfoil and boundary layers plot
% b_plotConv:   Show boundary layer convergence plot
function [d_star_t,d_star_b,I_crit,Cp,Cf,strengths,X_mid,Y_mid,norms,extra] = Thwaites_panel_1(X,Y,Vinf,alpha,mu,rho,it,b_plotFoil,b_plotConv)

nu = mu/rho;    % Kinematic Viscosity [m^2/s]
M = length(X);  % Number of points
N = M-1;        % Number of panels

% Set the geometry coordinates
X_base = X;
Y_base = Y;

% Instead of physically displacing the geometry due to the BL each iteration, a value is
% added to the freestream condition to simulate it
extra = zeros(N,1);

% Index of the separation point on the top, bottom, and overall
I_crit = zeros(3,1);

% Get the midpoints and normals from the geometry
[X_mid,Y_mid,norms,~,node_norms] = geo_decomp(X,Y);

% Perform an inviscid analysis to get the stagnation point
[Cp,strengths] = vortex_2order(X,Y,Vinf,alpha,extra,2*N);
[~,I] = max(Cp.*(X_mid<0.5));
I = I+1;

% Get the node normals and split them into top and bottom surfaces
%   (nodes are the endpoints of the panels)
norms_top_node = node_norms(I:end,:);
norms_bot_node = flip(node_norms(1:I,:));

BL_heights = zeros(it,1);

% Distance from the mids to sample the velocity
sample_dist = 0e-10;

% Preallocate panel skin friction coefficient vector 
Cf = zeros(N,2);

% 1. Panel method gives u at body
% 2. Solve for delta_star at those u
% 3. Get new body with those delta_star
% 4. Repeat 1-3
for k = 1:it
    if k ~= 1
        %% Inviscid Panel Method
        [Cp,strengths] = vortex_2order(X,Y,Vinf,alpha,extra,2*N);
        [X_mid,Y_mid,norms] = geo_decomp(X,Y);
    end
    
    % Get the normals at the panel mids for the top and bottom surfaces
    norms_top = norms(I:end,:);
    norms_bot = flip(norms(1:I,:));
    
    %% Create and analyze boundary layer on top & bottom surfaces
    
    % Split panel mids into top and bottom sections
    X_top = X_mid(I:end)';
    Y_top = Y_mid(I:end)';
    X_bot = flip(X_mid(1:I-1)');
    Y_bot = flip(Y_mid(1:I-1)');
    % Split nodes into top and bottom sections
    X_top_node = X(I:end)';
    Y_top_node = Y(I:end)';
    X_bot_node = flip(X(1:I)');
    Y_bot_node = flip(Y(1:I)');
    % Surface lengths
    S_top = 0;                                  % Total panel length so far (top surface)
    S_bot = 0;                                  % Total panel length so far (bottom surface)
    % Boundary layer heights
    deltas_top = zeros(length(X_top),2);        % Boundary layer height (top)
    delta_d_top = deltas_top;                   % Displacement boundary layer height (top)
    deltas_bot = zeros(length(X_bot),2);        % Boundary layer height (bottom)
    delta_d_bot = deltas_bot;                   % Displacement boundary layer height (bottom)
    
    Ue = zeros(length(X_top),1);    % Vector of surface velocities (AKA Ue, depending on the text)
    pos = Ue;                       % Vector of sampled positions along the surface (not x/c)
    lens = Ue;
    Re_t = Ue;                      % Vector of Reynolds Number (wrt momentum thickness) along surface
    
    %% Top surface
    % Getting surface velocity and location for every top panel
    for i = 1:length(X_top)
        % Control point
        xi = X_top(i);
        yi = Y_top(i);
        % Make variables for panel nodes for easier access
        dx = X_top_node(i + 1) - X_top_node(i);
        dy = Y_top_node(i + 1) - Y_top_node(i);
        len = sqrt(dx^2 + dy^2);    % Panel length
        lens(i) = len;
        
        % Store the position (length along surface)
        pos(i) = S_top + len/2;
        
        % getVel gets velocity at specified position in flowfield based on
        % inviscid solution
        Ue(i) = getVel(xi + norms_top(i,1)*sample_dist,yi + norms_top(i,2)*sample_dist,X,Y,strengths,Vinf,alpha);
        
        % Add to total length travelled
        S_top = S_top + len;
    end
    
    % Store the values for the last point
    pos(end) = S_top;
    
    % Calculate the velocity derivative from the velocity vector
    d_Ue = diff(Ue)./diff(pos);
    d_Ue = [d_Ue; d_Ue(end)];
    
    % Get the momentum thickness (theta) from the velocities
    %   theta^2 = 0.441*nu*Ue^-6 * int(Ue^5)
    theta = zeros(length(pos),1);
    for i = 1:length(Ue)
        if i == 1
            theta(i) = 0;
        else
            theta(i) = (0.441*nu*(Ue(i))^-6)*trapz(pos(1:i),Ue(1:i).^5);
            theta(i) = sqrt(theta(i));
        end
        % Momentum thickness Reynold Number
        Re_t(i) = Ue(i)*theta(i)/nu;
    end
    
    % Surface length Reynold Number
    Re_x = Ue.*pos./nu;
    
    % Get the initial value for the momentum thickness, theta0
    vd0 = d_Ue(1);
    vd1 = d_Ue(2);
    vdd0 = (vd1-vd0)/(pos(2)-pos(1));   % Second derivative of Ue (at stagnation point)
    theta(1) = theta(2)^2 - (-.45*mu/rho/7*vdd0/vd0^2)*(pos(2)-pos(1));
    theta(1) = sqrt(abs(theta(1)));
    theta(1) = 0;
    theta0 = theta(1);
    
    % Now the rest of the Thwaites-Pohlhausen method
    lambda = (theta.^2).*d_Ue/nu;       % Thwaites parameter
    % Pohlhausen functions g and f
    g_l = (lambda + 0.09).^0.62;
    f_l = 2 + 4.14*(0.25-lambda).^1 - 83.5*(0.25-lambda).^2 ...
        + 854*(0.25-lambda).^3 - 3337*(0.25-lambda).^4 ...
        + 4586*(0.25-lambda).^5;
    delta_d_star = f_l.*theta;              % Displacement thickness
    t_w_top = mu*Ue.*g_l./theta;  % Wall stress
    
    % Remove imaginary stresses (set to 0 if imaginary)
    for i=1:length(t_w_top)
        if ~isreal(t_w_top(i))
            t_w_top(i) = 0;
        end
    end
    
    Cft = t_w_top*2/rho/Vinf^2;     % Skin friction coefficient (top)
    Cft(isnan(Cft)) = 0;
    
    % Get separation point (where skin friction is 0 AKA when lambda is -0.09)
    x_crit = X_top_node(1:end-1).*(pos<interp1(lambda(1:end-1),pos(1:end-1),-.09));
    y_crit = Y_top_node(1:end-1).*(pos<interp1(lambda(1:end-1),pos(1:end-1),-.09));
    x_crit = x_crit(x_crit~=0);
    y_crit = y_crit(y_crit~=0);
    if isempty(x_crit)
        x_crit_top = X_top_node(end-1);
        y_crit_top = Y_top_node(end-1);
    else
        x_crit_top = x_crit(end);
        y_crit_top = y_crit(end);
    end
    % Separation point index in terms of top nodes array
    I_crit(1) = find(Y_top_node==y_crit_top,1);
    I_crit(3) = I_crit(1) + I;
    
    % Approximate the separated flow
    slope = (delta_d_star(I_crit(1))-delta_d_star(I_crit(1)-1))/(pos(I_crit(1))-pos(I_crit(1)-1));
    slope = 4*abs(slope);
    dds0 = delta_d_star(I_crit(1));
    for i=I_crit(1)+1:length(delta_d_star)
        delta_d_star(i) = dds0 + slope*(pos(i) - pos(I_crit(1)));
    end
    
    % Get the momentum thickness global coordinates for the top surface
    theta_top = zeros(length(theta),2);
    for i = 1:length(delta_d_top)
        theta_top(i,:) = theta(i)*norms_top_node(i,:) + [X_top(i) Y_top(i)];
    end
    
    % Convert displacement thickness from panel frames to fixed frame
    for i = 1:length(delta_d_top)
%         delta_d_top(i,:) = delta_d_star(i)*norms_top_node(i,:) + [X_top_node(i) Y_top_node(i)];
        delta_d_top(i,:) = delta_d_star(i)*norms_top_node(i,:) + [X_base(I+i-1)' Y_base(I+i-1)'];
    end
    
    % TE node
    delta_d_top = [delta_d_top; delta_d_star(length(delta_d_top))*norms_top_node(length(delta_d_top),:) + [X_base(I+length(delta_d_top))' Y_base(I+length(delta_d_top))']];
    
    d_star_t = delta_d_top;
    
    extra(I:end) = Ue.*([delta_d_star(2:end);delta_d_star(end)]-delta_d_star)./lens;
    
    % Clear for use on the bottom surface
    Ue = zeros(length(X_bot),1);
    lens = Ue;
    pos = Ue;
    
    %% Bottom surface
    % Getting surface velocity and location for every bottom panel
    for i = 1:length(X_bot)
        % Control point
        xi = X_bot(i);
        yi = Y_bot(i);
        % Make variables for panel nodes for easier access
        dx = X_bot_node(i + 1) - X_bot_node(i);
        dy = Y_bot_node(i + 1) - Y_bot_node(i);
        len = sqrt(dx^2 + dy^2);    % Panel length
        
        lens(i) = len;
        
        % Store the position (length along surface)
        pos(i) = S_bot + len/2;
        
        % getVel gets velocity at specified position in flowfield based on
        % inviscid solution
        Ue(i) = getVel(xi + norms_bot(i,1)*sample_dist,yi + norms_bot(i,2)*sample_dist,X,Y,strengths,Vinf,alpha);
        
        % Add to total length travelled
        S_bot = S_bot + len;
    end
    
    % Store the values for the last point
    pos(end) = S_bot;
    
    % Calculate the velocity derivative from the velocity vector
    d_Ue = diff(Ue)./diff(pos);
    d_Ue = [d_Ue; d_Ue(end)];
    
    % Get the momentum thickness, theta, from the velocities
    %   theta^2 = 0.441*nu*Ue^-6 * int(Ue^5)
    theta = zeros(length(pos),1);
    for i = 1:length(Ue)
        if i == 1
            theta(i) = 0;
        else
            theta(i) = (0.441*nu*(Ue(i))^-6)*trapz(pos(1:i),Ue(1:i).^5);
            theta(i) = sqrt(theta(i));
        end
    end
    theta(1) = theta0;
    
    % Now the rest of the Thwaites-Pohlhausen method
    lambda = (theta.^2).*d_Ue/nu;       % Thwaites parameter
    % Pohlhausen functions g and f
    g_l = (lambda + 0.09).^0.62;
    f_l = 2 + 4.14*(0.25-lambda).^1 - 83.5*(0.25-lambda).^2 ...
        + 854*(0.25-lambda).^3 - 3337*(0.25-lambda).^4 ...
        + 4586*(0.25-lambda).^5;
    delta_d_star = f_l.*theta;              % Displacement thickness
    t_w_bot = mu*Ue.*g_l./theta;  % Wall stress
    
    for i=1:length(t_w_bot)
        if ~isreal(t_w_bot(i))
            t_w_bot(i) = 0;
        end
    end
    
    Cfb = t_w_bot*2/rho/Vinf^2;     % Skin friction coefficient (bottom)
    Cfb(1) = 0;
    Cfb(isnan(Cfb)) = 0;
    Cf = {Cfb; Cft};        % Skin friction coefficient
    
    % Get separation point (where skin friction is 0 AKA when lambda is -0.09)
    x_crit = X_bot_node(1:end-1).*(pos<interp1(lambda(1:end-1),pos(1:end-1),-.09));
    y_crit = Y_bot_node(1:end-1).*(pos<interp1(lambda(1:end-1),pos(1:end-1),-.09));
    x_crit = x_crit(x_crit~=0);
    y_crit = y_crit(y_crit~=0);
    if isempty(x_crit)
        x_crit_bot = X_bot_node(end-1);
        y_crit_bot = Y_bot_node(end-1);
    else
        x_crit_bot = x_crit(end);
        y_crit_bot = y_crit(end);
    end
    % Separation point index in terms of bottom nodes array
    I_crit(2) = find(Y_bot_node==y_crit_bot,1);
    
    % Approximate the separated flow
    slope = (delta_d_star(I_crit(2))-delta_d_star(I_crit(2)-1))/(pos(I_crit(2))-pos(I_crit(2)-1));
    slope = abs(slope);
    dds0 = delta_d_star(I_crit(2));
    for i=I_crit(2)+1:length(delta_d_star)
        delta_d_star(i) = dds0 + slope*(pos(i) - pos(I_crit(2)));
    end
    
    % Get the momentum thickness global coordinates for the bottom surface
    theta_bot = zeros(length(theta),2);
    for i = 1:length(delta_d_bot)
        theta_bot(i,:) = theta(i)*norms_bot_node(i,:) + [X_bot(i) Y_bot(i)];
    end
    
    % Convert displacement thickness from panel frames to fixed frame
    for i = 1:length(delta_d_bot)
        delta_d_bot(i,:) = delta_d_star(i)*norms_bot_node(i,:) + [flip(X_base(I-i+1)') flip(Y_base(I-i+1)')];
    end
    
    % TE node
    delta_d_bot = [delta_d_bot; delta_d_star(length(delta_d_bot))*norms_bot_node(length(delta_d_bot),:) + [flip(X_base(1)') flip(Y_base(1)')]];
    
    d_star_b = delta_d_bot;
    
    extra(1:I-1) = flip(Ue.*([delta_d_star(2:end);delta_d_star(end)]-delta_d_star)./lens);
    
    col = 'b';
    if k==it
        col = 'm';
    end
    
    % Set new effective body to the displacement thickness
    X = [flip(delta_d_bot(:,1)); delta_d_top(2:end,1)]';
    Y = [flip(delta_d_bot(:,2)); delta_d_top(2:end,2)]';
    
end

% Get the final Cp and panel strengths from an inviscid method
[Cp,strengths,X_mid,Y_mid,norms] = vortex_2order(X,Y,Vinf,alpha,extra,I_crit(1)+I);

% Plot final boundary layer thicknesses and the stagnation point
if b_plotFoil
    figure(2)
    p(1) = plot(X_base,Y_base,'k.-');
    hold on
    p(2) = plot(deltas_top(:,1),deltas_top(:,2),'r');
    p(3) = plot(delta_d_top(:,1),delta_d_top(:,2),'b');
    plot(deltas_bot(:,1),deltas_bot(:,2),'r')
    plot(delta_d_bot(:,1),delta_d_bot(:,2),'b')
    p(4) = plot(X_base(I),Y_base(I),'m*');
    pbaspect([1.4 1.0 1])
    xlim([-.2 1.2])
    ylim([-.5 .5])
    xlabel("x/c")
    ylabel("y/c")
    title("Airfoil Boundary Layers")
    hold off
    legend([p(1) p(2) p(3) p(4)],'Geometry','\delta','\delta^*','Stagnation Point')
end

% Plot of displacement BL height at a station vs iteration number
if b_plotConv
    figure()
    plot(BL_heights)
    ylim([0 max(BL_heights)])
    xlabel("Iteration #")
    ylabel("\delta^*")
    title("\delta^* Convergence")
    grid on
end

end



% ==========================================================================
%               AUXILIARY FUNCTIONS
% ==========================================================================

%% Get induced velocity from vortices at a point (x,y) due to a panel with end points (xj,yj) and (xjj,yjj), and strengths gj and gjj
function [u1,v1,u2,v2] = getVi(gj,gjj,x,y,xj,yj,xjj,yjj,same)
a = atan2(yjj-yj,xjj-xj);   % Panel orientation angle [rad]

% Transformation Matrix
A = [cos(a) sin(a);
    -sin(a) cos(a)];

% Convert to panel coordinates
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

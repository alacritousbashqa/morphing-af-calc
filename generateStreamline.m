function [streamline] = generateStreamline(start,X,Y,strengths,Vinf,alpha,Nw)
%% Creates a streamline (matrix of points) from a starting point based on the velocity induced by the airfoil vortices and the freestream

N = Nw;                     % Number of wake panels
streamline = zeros(N,2);    % Wake from the separation point

% Panel lengths increase based on cosine distribution
t = linspace(0,pi,N);
initialPanelLength = 2e-5;  % Length of first wake panel
L = initialPanelLength + 0.05*(1 - cos(t));

streamline(1,:) = [start(1) start(2)];   % Set the first wake point at the separation point

% The panels should follow the streamlines of the air flow
dir_new = 0;
for i = 1:N-1
    dist = L(i)/2;    % Collocation point distance from previous wake node
    coll = [streamline(i,1)+dist*cos(dir_new) streamline(i,2)+dist*sin(dir_new)];   % Collocation point (AKA panel center)
    dir_old = atan2(coll(2)-streamline(i,2),coll(1)-streamline(i,1));   % Initial panel angle
    % Iterate until it follows streamline (within 1e-4 angle difference between iterations) or 100 iterations are done
    for j = 1:100
        [~,u,v] = getVel(coll(1),coll(2),X,Y,strengths,Vinf,alpha); % Get velocity at collocation point
        dir_new = atan2(v,u);   % Get new direction based on new velocity calc
        if abs(dir_new-dir_old) <= 1e-4
            streamline(i+1,:) = [streamline(i,1)+dist*cos(dir_new) streamline(i,2)+dist*sin(dir_new)];
            break
        else
            coll = [streamline(i,1)+dist*cos(dir_new) streamline(i,2)+dist*sin(dir_new)];
            dir_old = atan2(coll(2)-streamline(i,2),coll(1)-streamline(i,1));
        end
    end
end

end
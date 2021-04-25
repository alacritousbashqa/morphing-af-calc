function [streamline] = generateStreamline(start,bCosine,pLength,X,Y,strengths,Vinf,alpha,Nw)
%% Creates a streamline (matrix of points) from a starting point based on the velocity induced by the airfoil vortices and the freestream

streamline = zeros(Nw,2);    % Wake from the separation point

if bCosine
    % Panel lengths increase based on cosine distribution
    t = linspace(0,pi,Nw);
    L = pLength + 0.05*(1 - cos(t));
else
    % Panel lengths remain constant
    L = pLength.*ones(Nw,1);
end

streamline(1,:) = [start(1) start(2)];   % Set the first wake point at the separation point

% The panels should follow the streamlines of the air flow
dir_new = 0;
for i = 1:Nw-1
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
        elseif j == 100
            streamline(i+1,:) = [streamline(i,1)+dist*cos(dir_old) streamline(i,2)+dist*sin(dir_old)];
        else
            coll = [streamline(i,1)+dist*cos(dir_new) streamline(i,2)+dist*sin(dir_new)];
            dir_old = atan2(coll(2)-streamline(i,2),coll(1)-streamline(i,1));
        end
    end
end

end
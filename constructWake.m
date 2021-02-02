function [spWake,teWake] = constructWake(SP,TE,X,Y,strengths,Vinf,alpha)

N = 500;                % Number of wake panels
spWake = zeros(N,2);    % Wake from the separation point
teWake = spWake;        % Wake from the trailing edge

% Panel lengths increase based on cosine distribution
t = linspace(0,pi,N);
L = 0.05*(1 - cos(t));

dist = 1e-5;    % Wake starting distance from surface point
spWake(1,:) = [SP(1)+dist SP(2)];   % Set the first wake just next to the separation point
teWake(1,:) = [TE(1)+dist TE(2)];   % "             "               "     trailing edge

% The panels should follow the streamlines of the air flow
for i = 1:N
    % Get the velocity
    [V,u,v] = getVel(spWake(i,1),spWake(i,1),X,Y,strengths,Vinf,alpha);
    % Get the scaled direction
    dir = L(i)*[u v]/V;
    spWake(i+1,:) = spWake(i,:) + dir;
    
    [V,u,v] = getVel(teWake(i,1),teWake(i,1),X,Y,strengths,Vinf,alpha);
    dir = L(i)*[u v]/V;
    teWake(i+1,:) = teWake(i,:) + dir;
end

end
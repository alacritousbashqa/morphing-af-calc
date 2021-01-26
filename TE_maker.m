%% Generate different trailing edge based on new point given in input TE
%   bot_start: index to start new TE on bottom surface
%   top_start: index to start new TE on top surface
%   TE: new trailing edge point
%   N_t: number of nodes to use in new top surface
%   N_b: number of nodes to use in new bottom surface
function [bot,top] = TE_maker(bot_start,top_start,TE,N_t,N_b)

% Generate linear function to calculate points
y_t = @(x) ((TE(2)-top_start(2))/(TE(1)-top_start(1)))*(x - top_start(1)) + top_start(2);
y_b = @(x) ((TE(2)-bot_start(2))/(TE(1)-bot_start(1)))*(x - bot_start(1)) + bot_start(2);

% Generate x points
x_t = linspace(top_start(1),TE(1),N_t);
x_b = linspace(TE(1),bot_start(1),N_b);

% Calculate (x,y)
top = [x_t; y_t(x_t)];
bot = [x_b; y_b(x_b)];

end
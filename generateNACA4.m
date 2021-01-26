%% Generate 4 digit NACA coordinates for n panels
function [X,Y] = generateNACA4(digits, n)
% Generate coordinates if it is a circle
if digits == "circle"
    t = linspace(0,2*pi,n+1)';
    X = cos(t);
    Y = sin(t);
    return;
end

% Get parameters from digits, assuming it was four digits
t = str2double(extractBetween(digits,3,4))/100;
pc = str2double(extractBetween(digits,2,2))/10;
mc = str2double(extractBetween(digits,1,1))/100;

if mod(n,2) == 1
    n = n/2+.5;
else
    n = round((n+1)/2);
end
beta = linspace(0,pi,n);
x = (1-cos(beta))/2;
x1 = 0;
x2 = 0;
for i=1:length(x)
    if x(i)<pc
        x1(i) = x(i);
    else
        x2(i) = x(i);
    end
end
x2(1:length(x1)) = [];

% y and dy/dx of camberline
yc = horzcat(mc/pc^2.*(2*pc.*x1-x1.^2),mc/(1-pc)^2.*(1-2*pc+2*pc.*x2-x2.^2));
dyc = horzcat(2*mc/pc^2.*(pc-x1),2*mc/(1-pc)^2.*(pc-x2));

% Thickness
a0 = .2969;
a1 = -.126;
a2 = -.3537;
a3 = .2843;
a4 = -.1015;
yt = t/0.2*(a0.*x.^(.5)+a1.*x+a2.*x.^2+a3.*x.^3+a4.*x.^4);

theta = atan(dyc);

xu = x-yt.*sin(theta);
xl = x+yt.*sin(theta);
xl(1)=[];
yu = yc+yt.*cos(theta);
yl = yc-yt.*cos(theta);
yl(1)=[];

X = horzcat(flip(xu),xl);
Y = horzcat(flip(yu),yl);

for i=1:length(X)
    if isnan(X(i))
        X(i) = 0;
    end
    if isnan(Y(i))
        Y(i) = 0;
    end
end
end
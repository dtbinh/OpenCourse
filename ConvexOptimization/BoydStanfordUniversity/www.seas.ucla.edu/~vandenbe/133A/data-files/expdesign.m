function [A1, A2, A3] = expdesign;

% experiment design

D = 4;

% source positions
x1=[-1.5; 0];  x2=[-0.5; 0]; x3 = [0.5; 0]; x4 = [1.5; 0];

% sensors positions for 1st configuration
y1=[-1.5; D];  y2=[-0.5; D]; y3 = [0.5; D]; y4 = [1.5; D];
A1 = [1/norm(x1-y1)^2  1/norm(x2-y1)^2 1/norm(x3-y1)^2 1/norm(x4-y1)^2;
    1/norm(x1-y2)^2  1/norm(x2-y2)^2 1/norm(x3-y2)^2 1/norm(x4-y2)^2;
    1/norm(x1-y3)^2  1/norm(x2-y3)^2 1/norm(x3-y3)^2 1/norm(x4-y3)^2;
    1/norm(x1-y4)^2  1/norm(x2-y4)^2 1/norm(x3-y4)^2 1/norm(x4-y4)^2];

% sensor positions for 2nd configuration
y1=[-1.5; D];  y2=[1.5; D]; y3 = [-1.5; -D]; y4 = [1.5; -D];
A2 = [1/norm(x1-y1)^2  1/norm(x2-y1)^2 1/norm(x3-y1)^2 1/norm(x4-y1)^2;
    1/norm(x1-y2)^2  1/norm(x2-y2)^2 1/norm(x3-y2)^2 1/norm(x4-y2)^2;
    1/norm(x1-y3)^2  1/norm(x2-y3)^2 1/norm(x3-y3)^2 1/norm(x4-y3)^2;
    1/norm(x1-y4)^2  1/norm(x2-y4)^2 1/norm(x3-y4)^2 1/norm(x4-y4)^2];

% sensor positions for 3rd configuration
y1=[-1.5; D];  y2=[1.5; D]; y3 = [-1.5-D; 0]; y4 = [1.5+D; 0];
A3 = [1/norm(x1-y1)^2  1/norm(x2-y1)^2 1/norm(x3-y1)^2 1/norm(x4-y1)^2;
    1/norm(x1-y2)^2  1/norm(x2-y2)^2 1/norm(x3-y2)^2 1/norm(x4-y2)^2;
    1/norm(x1-y3)^2  1/norm(x2-y3)^2 1/norm(x3-y3)^2 1/norm(x4-y3)^2;
    1/norm(x1-y4)^2  1/norm(x2-y4)^2 1/norm(x3-y4)^2 1/norm(x4-y4)^2];


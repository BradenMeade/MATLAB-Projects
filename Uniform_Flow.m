% Clears all information from previous runs
clear variables
close all

% Define all parameters

U = 1;                  % Free stream velocity

Alpha = 0*(pi/180);     % Flow angle of free stream

m = 30;                 % Length of iteration

% Defines x and y directions as 1xm vectors by specifying bounds and increment

x = linspace(-10,10,m);

y = linspace(-10,10,m);

% Meshgrid function generates mxm matrix containing all coordinate pairs obtained from 1xm vectors above 

[X,Y] = meshgrid(x,y);

% Scalar potential function which produces an mxm matrix

phi = U*(X*cos(Alpha)+Y*sin(Alpha));

% Velocity components of uniform flow. Note: Here, the 'ones' function is
% incorporated to return a constant-valued mxm matrix

v_x = U*cos(Alpha)*ones(m);

v_y = U*sin(Alpha)*ones(m);

% Generates uniform flow scalar potential field
figure;
contour(x,y,phi,'ShowText','on')
title('Uniform Flow Scalar Potential Field')

% Generates uniform flow velocity vector field
figure;
quiver(X,Y,v_x,v_y,0.5);
title('Uniform Flow Velocity Vector Field')
xlim([-5,5])
ylim([-5,5])
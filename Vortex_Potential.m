% Clears and closes any information from previous runs
clear variables
close all

% Define all mathematical and physical parameters to be used

R = 1;                              % Radius of cylinder

Gamma = 12;                         % Circularity of flow

m = 30;                             % Increment of iteration

r = linspace(R, 5, m);              % Radial distance respective to center of cylinder

Theta = linspace(0,2*pi, m);        % Rotation from parallel of free stream velocity direction


% Pre-allocating dependant variables to blank mxm matricies (program performance only)

U_x = zeros(m,m);

U_y = zeros(m,m);

x = zeros(m,m);

y = zeros(m,m);

U_r = zeros(m,m);

U_t = zeros(m,m);

Pot_function = zeros(m,m);

Stream_function = zeros(m,m);


for i = 1:m

    for j = 1:m

        % Basic definitions for Cartesian representation of polar coordinates

        x(i,j) = r(i)*cos(Theta(j));
        
        y(i,j) = r(i)*sin(Theta(j));

        % Potential function for source/sink

        Pot_function(i,j) = (Gamma/(2*pi))*Theta(j);

        % Generates the circumferential velocity component

        U_t(i,j) = (1/r(i))*(Gamma/(2*pi));

        U_r = zeros(m,m);

        % Convert radial to Cartesian velocity components

        U_x(i,j) = -sin(Theta(j))*U_t(i,j)+cos(Theta(j))*U_r(i,j);

        U_y(i,j) = sin(Theta(j))*U_r(i,j)+cos(Theta(j))*U_t(i,j);

    end

end


% Generates
figure;
contour(x,y,Pot_function,0.5*m, 'ShowText','on')
hold on
cylinder(R)
hold off
axis equal
axis tight
xlim([-3,3])
ylim([-3,3])
title('Irrotational Vortex Scalar Potential Field')

figure;
quiver(x,y,U_x, U_y, 0.5)
hold on
cylinder(R)
hold off
axis equal
axis tight
xlim([-3,3])
ylim([-3,3])
title('Irrotational Vortex Velocity Field')
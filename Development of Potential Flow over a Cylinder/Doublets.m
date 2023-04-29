% Clears and closes any information from previous runs
clear variables
close all

% Define all mathematical and physical parameters to be used
R = 1;

m = 200;                             % Increment of iteration

r = linspace(1, 10, m);              % Radial distance respective to center of cylinder

Theta = linspace(0,2*pi, m);        % Rotation from parallel of free stream velocity direction

K = 10;                             % Magnitude of doublet

% Pre-allocating dependant variables to blank mxm matricies (program performance only)

x = zeros(m,m);

y = zeros(m,m);

Pot_function = zeros(m,m);

Stream_function = zeros(m,m);


for i = 1:m

    for j = 1:m

        % Basic definitions for Cartesian representation of polar coordinates

        x(i,j) = r(i)*cos(Theta(j));
        
        y(i,j) = r(i)*sin(Theta(j));

        % Potential function for doublet

        Pot_function(i,j) = (K*cos(Theta(j)))/r(i);

        % Stream function for doublet

        Stream_function(i,j) = (-K*sin(Theta(j)))/r(i);

    end

end
v = [-7,-6,-5,-4,-3,0,3,4,5,6,7];

% Generates contour plot of potential function doublet
figure;
contour(x,y,Pot_function, v, 'ShowText','on')
hold on
cylinder(R)
hold off
axis equal
axis tight
xlim([-3.5,3.5])
ylim([-3.5,3.5])
title('Doublet Scalar Potential Function')

% Generates contour plot of stream function doublet
figure;
contour(x,y,Stream_function, v, 'ShowText','on')
hold on
cylinder(R)
hold off
axis equal
axis tight
xlim([-3.5,3.5])
ylim([-3.5,3.5])
title('Doublet Stream Function')

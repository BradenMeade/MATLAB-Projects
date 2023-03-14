% Clears and closes any information from previous runs
clear variables
close all

% Define all mathematical and physical parameters to be used

R = 1;                              % Radius of cylinder

m = 30;                             % Increment of iteration

r = linspace(0.5, 15, m);           % Radial distance respective to center of cylinder

Theta = linspace(0,2*pi, m);        % Rotation from parallel of free stream velocity direction

Q = 1;                              % Strength of source/sink (volumetric flow rate)


% Pre-allocating variables to blank mxm matricies (program performance only)

U_x = zeros(m,m);

U_y = zeros(m,m);

x = zeros(m,m);

y = zeros(m,m);

U_r = zeros(m,m);

U_t = zeros(m,m);

Pot_function = zeros(m,m);

Stream_function = zeros(m,m);


% For loop to compute dependant variables with iteration from range 1:m

for i = 1:m

    for j = 1:m

        % Basic definitions for Cartesian representation of polar coordinates

        x(i,j) = r(i)*cos(Theta(j));
        
        y(i,j) = r(i)*sin(Theta(j));

        % Potential function for source/sink

        Pot_function(i,j) = (Q/(2*pi))*log(r(i)/R);

        % Radial and tangential velocity components

        U_r(i,j) = Q/(2*pi*r(i));

        U_t = zeros(m,m);

        % Convert radial to Cartesian velocity components

        U_x(i,j) = -sin(Theta(j))*U_t(i,j)+cos(Theta(j))*U_r(i,j);

        U_y(i,j) = sin(Theta(j))*U_r(i,j)+cos(Theta(j))*U_t(i,j);

    end

end

% Specifies contour elevations to be plotted

v = [-0.15,-0.1,0,0.05,0.1,0.15,0.2];

% Generates source potential contour plot
figure;
contour(x,y,Pot_function,v,'ShowText','on')
hold on
plot([-5,5],[1,1],"--r")
plot([-5,5],[-1,-1],"--r")
hold off
axis equal
axis tight
xlim([-3,3])
ylim([-3,3])
title('Source/Sink Scalar Potential Field')

% Generates source velocity vector plot
figure;
quiver(x,y,U_x, U_y,0,v)
hold on
cylinder(R)
hold off
axis equal
axis tight
xlim([-3,3])
ylim([-3,3])
title('Source/Sink Velocity Vector Field')

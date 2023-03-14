% Clears and closes any information from previous runs
clear variables
close all

% Define all mathematical and physical parameters to be used

R = 1;                              % Radius of cylinder

U = 1;                              % Free stream velocity (m/s)

P = 101000;                         % Free stream pressure (Pa)

rho = 997;                          % Fluid density (kg/m^3)

Gamma = -12;                        % Circularity of flow

m = 30;                             % Increment of iteration

r = linspace(R, 5*R, m);            % Radial distance respective to center of cylinder

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

for i=1:m

    for j=1:m

        % Basic definitions for Cartesian representation of polar coordinates

        x(i,j) = r(i)*cos(Theta(j));
        
        y(i,j) = r(i)*sin(Theta(j));
        
        % Radial and tangential velocity components of rotating cylinder

        U_r(i,j) = U*(1-(R/r(i))^2)*cos(Theta(j));

        U_t(i,j) = -U*(1+(R/r(i))^2)*sin(Theta(j))+(Gamma/(2*pi*r(i)));

        % Cartesian conversion

        U_x(i,j) = -sin(Theta(j))*U_t(i,j)+cos(Theta(j))*U_r(i,j);

        U_y(i,j) = sin(Theta(j))*U_r(i,j)+cos(Theta(j))*U_t(i,j);

        % Velocity magnitude

        U_abs = sqrt((U_x).^2+(U_y).^2);

        % Definition of stream function for rotating cylinder

        Stream_function(i,j) = U*sin(Theta(j))*(r(i)-(R^2/r(i)))-(Gamma/(2*pi))*log(r(i)/R);

        % Pressure gradient

        p = (1/2)*rho*(U.^2-(U_abs).^2)+P;

    end

end


% Generates stream function plot
figure;
contour(x,y,Stream_function, m);
hold on
cylinder(R)
hold off
axis equal
grid on
title 'Stream Function Contour - Rotating Cylinder'
xlim([-3,3])
ylim([-3,3])

% Generates velocity field
figure;
quiver(x,y,U_x,U_y);
hold on
cylinder(R)
hold off
axis equal
title('Velocity Vector Field - Rotating Cylinder')
xlim([-3,3])
ylim([-3,3])

% Generates velocity graident plot
figure;
vmap = pcolor(x,y,U_abs);
set(vmap,'EdgeColor','none')
axis equal
grid off
xlim([-3,3])
ylim([-3,3])
colorbar
shading interp
title 'Velocity Gradient - Rotating Cylinder'

% Generates pressure gradient plot
figure;
pmap = pcolor(x,y,p);
set(pmap, 'EdgeColor','none')
axis equal
xlim([-3,3])
ylim([-3,3])
colorbar
shading interp
title 'Pressure Gradient - Rotating Cylinder'

% Generates layered stream function & velocity vector plot
figure;
contour(x,y,Stream_function, m);
hold on
quiver(x,y,U_x, U_y)
hold on
cylinder(R)
hold off
axis equal
axis tight
xlim([-3,3])
ylim([-3,3])

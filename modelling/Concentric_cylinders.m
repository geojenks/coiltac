%% 2D Initialisation

% Units correspond to 1/10 mm
% Enter the dimensions
Nx = 501;     % Number of X-grids
Ny = 501;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
   
Ni = 10000;  % Number of iterations for the Poisson solver
% comment out to use input from last time
V = zeros(Nx,Ny);   % Potential (Voltage) matrix
error = [];
T = 0;            % Top-wall potential
B = 0;            % Bottom-wall potential
L = 0;            % Left-wall potential
R = 0;            % Right-wall potential
%-------------------------------------------------------------------------%
% Initializing edges potentials
%-------------------------------------------------------------------------%

V(1,:) = L;
V(Nx,:) = R;
V(:,1) = B;
V(:,Ny) = T;
%-------------------------------------------------------------------------%
% Initializing Corner potentials
%-------------------------------------------------------------------------%
V(1,1) = 0.5*(V(1,2)+V(2,1));
V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));
V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));

%% 2D
%% wire cross sections (concentric wires)

wire_1_radius = 100;
wire_2_radius = 200;
wire_1_centre = [mpx, mpy];
wire_2_centre = [mpx, mpy];
err_start = length(error);
check_every = 10;
error = [error zeros(1,Ni/check_every)];
perc = 0;
for z = 1:Ni    % Number of iterations
    if round(100*z/Ni) ~= perc
        disp(num2str(round(100*z/Ni)) + "%")
        perc = round(100*z/Ni);
    end
    for i=1:Nx
        for j=1:Ny
            % The next two lines are to force the matrix to hold the 
            % potential values for all iterations
            if (i-(wire_1_centre(1)))^2 + (j-(wire_1_centre(2)))^2 < wire_1_radius^2
                V(i,j) = 100;
            elseif (i-(wire_2_centre(1)))^2 + (j-(wire_2_centre(2)))^2 > wire_2_radius^2
                V(i,j) = -100;
            else
                V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
            end
        end
    end
    if mod(z,check_every) == 0
        % Take transpose for proper x-y orientation
        V = V';
        [Ex,Ey]=gradient(V);
        Ex = -Ex;
        Ey = -Ey;
        % Electric field Magnitude
        E = sqrt(Ex.^2+Ey.^2);
        perm = .00000000000885;
        Volt=200;
        Field = mean(nonzeros(elements_within_radius(E,[mpx, mpy],[wire_1_radius+1 wire_1_radius+2])));
        %figure;plot(nonzeros(results))
        %mean(nonzeros(results))
        density = perm * Field;
        Q = density * 2*pi*wire_1_radius;
        numericalCap = Q/Volt;
        %analyticalCap = perm*pi/(acosh((2*position_wire)/(2*wire_1_radius)));
        analyticalCap = 2*pi*perm/ log(wire_2_radius/wire_1_radius);
        error(err_start+z/check_every) = analyticalCap - numericalCap;
        V=V';
    end
end


%% 2D plotting etc.:
% Take transpose for proper x-y orientation

V = V';
[Ex,Ey]=gradient(V);
Ex = -Ex;
Ey = -Ey;
% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);  
x = (1:Nx)-mpx;
y = (1:Ny)-mpy;
% Contour Display for electric potential

figure(1)
contour_range_V = -101:0.5:101;
contour(x,y,V,contour_range_V,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
title('Electric Potential distribution, V(x,y) in volts','fontsize',14);
h1=gca;
set(h1,'fontsize',14);
fh1 = figure(1); 
set(fh1, 'color', 'white')
% Contour Display for electric field
figure(2)
contour_range_E = -20:0.05:20;
contour(x,y,E,contour_range_E,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
title('Electric field distribution, E (x,y) in V/m','fontsize',14);
h2=gca;
set(h2,'fontsize',14);
fh2 = figure(2); 
set(fh2, 'color', 'white')

% Quiver Display for electric field Lines
figure(3)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,Ex,Ey,2)%,'filled','.')
title('Electric field Lines, E (x,y) in V/m','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(3); 
set(fh3, 'color', 'white')
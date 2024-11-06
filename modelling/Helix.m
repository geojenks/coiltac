%% 3D Initialisation

perm = .00000000000885;
Volt = 200;

Nx = 61;     % Number of X-grids
Ny = 61;     % Number of Y-grids
Nz = 61;     % Number of Z-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
mpz = ceil(Nz/2); % Mid point of z

Ni = 500;  % Number of iterations for the Poisson solver

%% Cylindrical core and helix configuration
% Units are 1/10 mm

d = 15; %radial distance centrepoint of wire
w = 5; % radial thickness of wire
r = 5; % radial thickness of core
pitches = [60];
% if multiple pitches needed:
%pitches = [1, 10,20,60,100,140,180, 220, 300, 500];


for q = [1] % can iterate through e.g. [1:10] if varying pitch height
    pitch = pitches(q); %60
    % make the cylinder protrude 3/4 of the way into the helix
    cyl_height = 3*Nz/4;
    theta = atan2((2*pi*d),pitch); % static pitch angle (theta can be input directly)
    %for theta = 0:0.1:pi/2-0.01 % to plot the changing cross section
    % the eccentric radius is w+w*tan(theta), where w is the circular radius
    eccentric = w+w*sin(theta);
    % the angle that gives this as an arc length
    psimax = (w+w*tan(theta))/d;

    % Build ellipse with eccentricity defined by the pitch angle
    ellipse = []; % x and y coordinates of ellipse points
    samples = 50;
    for angle = pi/2:pi/(samples/2):3*pi/2
        ellipse = [ellipse; w*cos(angle),   (w+w*sin(theta))*sin(angle),...
                            w*cos(angle+pi),(w+w*sin(theta))*sin(angle+pi)];
        % we only really need 2 columns, and to dupplicate 1 with a minus
    end
    % append angles to map the heights to
    ellipse(:,5) = ellipse(:,2)/max(ellipse(:))*psimax;
    % set voltage
    V = zeros(Nx,Ny,Nz);
    perc = 0;

    for it = 1:Ni    % Number of iterations
        % Display progress every 5%
        if round(100*it/Ni) ~= perc
            if mod(perc, 5) == 0
                disp(num2str(round(100*it/Ni)) + "%")
            end
            perc = round(100*it/Ni);
        end
        % Laplacian smoothing
        for k=1:Nz
            turn = mod(pi+ k*2*pi/pitch, 2*pi) - pi; % turns calculated from z change and pitch
            for i=1:Nx % this is actually y
                for j=1:Ny %this is actually x
                    % The next two lines are to force the matrix to hold the 
                    % potential values for all iterations. Derichlet
                    % boundary conditions
                    if (sqrt((i-mpx)^2+(j-mpy)^2) - d)^2 + (((mod(atan2((j-mpy),(i-mpx))+ turn +pi, 2*pi)-pi)*ellipse(10,2)/ellipse(10,5))*w/eccentric)^2 < w^2
                        V(i,j,k) = Volt;
                    elseif (i-mpx)^2 + (j-mpy)^2 < r^2 && k < cyl_height
                        V(i,j,k) = 0;
                    else
                        valid_neighbors = 0;
                        neighbor_sum = 0;
                        for ii = -1:1
                            for jj = -1:1
                                for kk = -1:1
                                    % Do not include the central cell itself
                                    if (ii == 0) && (jj == 0) && (kk == 0)
                                        continue;
                                    end
                                    % Check if neighbor is within the valid range
                                    if (i+ii >= 1) && (i+ii <= Nx) && (j+jj >= 1) && (j+jj <= Ny) && (k+kk >= 1) && (k+kk <= Nz)
                                        valid_neighbors = valid_neighbors + 1;
                                        neighbor_sum = neighbor_sum + V(i+ii, j+jj, k+kk);
                                    end
                                end
                            end
                        end
                        % Calculate the average of the valid neighbors
                        V(i, j, k) = neighbor_sum / valid_neighbors;
                    end
                end
            end
        end
    end

    %For each "slice" work out E-field
    [i_grid, j_grid] = meshgrid(1:Nx, 1:Ny);
    squared_distances = (i_grid - mpy).^2 + (j_grid - mpx).^2;
    mask_core = squared_distances < (r+1)^2 & squared_distances >= r^2;
    % Initialise E-Field
    Exyx = zeros(Nx,Ny,Nz);
    Exyy = zeros(Nx,Ny,Nz);
    Ez = zeros(Nx,Ny,Nz);
    E = zeros(Nx,Ny,Nz);
    for k = 2:ceil(cyl_height)
        if k == ceil(cyl_height)
            % the top of the cylinder is a bit different
            mask_core = squared_distances < (r+1)^2;
        end
        [Exyx(:,:,k), Exyy(:,:,k)] = gradient(V(:,:,k)');
        Exyx(:,:,k) = -Exyx(:,:,k) .* mask_core;
        Exyy(:,:,k) = -Exyy(:,:,k) .* mask_core;
        Ez(:,:,k) = ((V(:,:,k-1)-V(:,:,k+1))/2) .* mask_core;
        E(:,:,k) = sqrt(Exyx(:,:,k).^2+Exyy(:,:,k).^2+Ez(:,:,k).^2);
        % Model field as homogeneously spread over surface
        Field(k) = mean(nonzeros(E(:,:,k)));
        % calculate this idealised density and charge
        density(q,k) = perm * Field(k);
        Q(k) = density(q,k) * pi * 2 * r * cyl_height;
        % estimate unit capacitance
        numericalCap(q,k) = Q(k)/Volt;
    end
    % convert to SI units
    Capacitance(q) = sum(numericalCap(q,:));
    conc_pred(q) = 2*pi*perm/ log((d-w)/r)*(cyl_height);%(mm)
    par_pred(q)  = perm*pi/(acosh(((d+2*r))/(r*2)))*(cyl_height);
end

%% 3D plotting
% we can only really plot a slice of the data, offset to see up or down
offset = -20;

% in-plane gradient
[Exyx, Exyy] = gradient(V(:,:,mpz+offset)');
Exyx = -Exyx;
Exyy = -Exyy;
Ez = (V(:,:,mpz+offset-1)-V(:,:,mpz+offset+1))/2;

% Electric field Magnitude
E = sqrt(Exyx.^2+Exyy.^2+Ez.^2);
x = (1:Nx);%-mpx;
y = (1:Ny);%-mpy;
Z = mpz+offset; % plot the central z-slice
% Contour Display for electric potential
figure(1)
contour_range_V = -101:0.5:101;
contour(V(:,:,Z),contour_range_V,'linewidth',0.5);
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
contour(E(:,:),'linewidth',0.5);
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
hold on,
quiver3(zeros([Nx, Ny]),Exyx, Exyy, Ez)

title('Electric field Lines, E (x,y) in V/m','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(3); 
set(fh3, 'color', 'white')

figure;
surf(Ez)
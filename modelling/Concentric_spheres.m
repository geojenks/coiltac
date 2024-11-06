%% 3D Initialisation

% Enter the dimensions
perm = .00000000000885;
Volt = 200;

Nx = 301;     % Number of X-grids
Ny = 301;     % Number of Y-grids
Nz = 301;     % Number of Z-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
mpz = ceil(Nz/2); % Mid point of z
   
Ni = 10000;  % Number of iterations for the Poisson solver

V = zeros(Nx,Ny,Nz);   % Potential (Voltage) matrix
error = [];
T = 0;            % Top-wall potential
B = 0;            % Bottom-wall potential
L = 0;            % Left-wall potential
R = 0;            % Right-wall potential
U = 0;            % Upper-wall potential
D = 0;            % Lower/Down-wall potential
%-------------------------------------------------------------------------%
% Initializing edges potentials
%-------------------------------------------------------------------------%
V(1,:,:) = L;
V(Nx,:,:) = R;
V(:,1,:) = B;
V(:,Ny,:) = T;
V(:,:,1) = D;
V(:,:,Nz) = U;
%-------------------------------------------------------------------------%
% Initializing Corner potentials
%-------------------------------------------------------------------------%
V(1,1,1) = 0.5*(V(1,2)+V(2,1));
V(Nx,1,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
V(1,Ny,1) = 0.5*(V(1,Ny-1)+V(2,Ny));
V(Nx,Ny,1) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));
V(1,1,Nz) = 0.5*(V(1,2)+V(2,1));
V(Nx,1,Nz) = 0.5*(V(Nx-1,1)+V(Nx,2));
V(1,Ny,Nz) = 0.5*(V(1,Ny-1)+V(2,Ny));
V(Nx,Ny,Nz) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));%}

%% Spherical 3D

core_radius = 75; % from 100
col = 3; %which column to put the data in
%outer_radius = 90;
%numericalCap = zeros([2,3]);
%analyticalCap = zeros([2,3]);
%errors = zeros([2,3]);
for outer_radius = 150%, 80, 90]
    err_start = length(error);
    %check_every = 10;
    %error = [error zeros(1,Ni/check_every)];
    perc = 0;
    Nx = 2* outer_radius + 3;
    Ny = 2* outer_radius + 3;
    Nz = 2* outer_radius + 3;
    mpx = ceil(Nx/2); % Mid-point of x
    mpy = ceil(Ny/2); % Mid point of y
    mpz = ceil(Nz/2); % Mid point of z
    core_centre = [mpx, mpy, mpz];
    V = zeros(Nx,Ny,Nz);
    % Create meshgrid for i, j, k indices
    [i_grid, j_grid, k_grid] = meshgrid(1:Nx, 1:Ny, 1:Nz);
    % Calculate squared distances from the core_centre
    squared_distances = (i_grid - core_centre(1)).^2 + (j_grid - core_centre(2)).^2 + (k_grid - core_centre(3)).^2;
    % Create masks for the conditions in the if-elseif statements
    mask_core = squared_distances < core_radius^2;
    mask_outer = squared_distances > outer_radius^2;
    % Set values for the core and outer regions
    V(mask_core) = 0;
    V(mask_outer) = Volt;
    % Create a mask for the remaining region
    mask_remaining = ~mask_core & ~mask_outer;
    % Find the indices of the masked elements
    [ind_i, ind_j, ind_k] = ind2sub(size(mask_remaining), find(mask_remaining));
    for z = 1:Ni    % Number of iterations
        if round(100*z/Ni) ~= perc
            if mod(perc, 5) == 0
                disp(num2str(round(100*z/Ni)) + "%")
            end
            perc = round(100*z/Ni);
        end
        % Iterate through the masked values
        for idx = 1:length(ind_i)
            i = ind_i(idx);
            j = ind_j(idx);
            k = ind_k(idx);
            V(i,j,k)=(1/6)*(V(i+1,j,k)+V(i-1,j,k)+V(i,j+1,k)+V(i,j-1,k)+V(i,j,k+1)+V(i,j,k-1)); %
        end
        % Calculate values for the remaining region
        % Create a temporary copy of V
        %V_temp = zeros(Nx, Ny, Nz);
        % Calculate values for the remaining region using V_temp
        %V_temp(2:(Nx-1), 2:(Ny-1), 2:(Nz-1)) = (1/6) * (V_temp(3:Nx, 2:(Ny-1), 2:(Nz-1)) + V_temp(1:(Nx-2), 2:(Ny-1), 2:(Nz-1)) + V_temp(2:(Nx-1), 3:Ny, 2:(Nz-1)) + V_temp(2:(Nx-1), 1:(Ny-2), 2:(Nz-1)) + V_temp(2:(Nx-1), 2:(Ny-1), 3:Nz) + V_temp(2:(Nx-1), 2:(Ny-1), 1:(Nz-2)));
        % Now update V with the calculated values from V_temp for the remaining region
        %V(~mask_core & ~mask_outer) = V_temp(~mask_core & ~mask_outer);
        %{
        for k=1:Nz
            for i=1:Nx
                for j=1:Ny   
                    % The next two lines are to force the matrix to hold the 
                    % potential values for all iterations
                    if (i-(core_centre(1)))^2 + (j-(core_centre(2)))^2 + (k-(core_centre(3)))^2 < core_radius^2
                        V(i,j,k) = 0;
                    elseif (i-(core_centre(1)))^2 + (j-(core_centre(2)))^2 + (k-(core_centre(3)))^2 > outer_radius^2
                        V(i,j,k) = Volt;
                    else
                        V(i,j,k)=(1/6)*(V(i+1,j,k)+V(i-1,j,k)+V(i,j+1,k)+V(i,j-1,k)+V(i,j,k+1)+V(i,j,k-1)); %
                    end
                end
            end
        end
        %}
        %{
        if mod(z,check_every) == 0
            E = zeros(Nx,Ny,Nz);
            % Can be optimised by making it only calculate E within the desired
            % radius
            for k = 2:Nz-1
                [Exyx, Exyy] = gradient(V(:,:,k)');%-gradient(permute(conj(V),[1,2,3]));
                Exyx = -Exyx;
                Exyy = -Exyy;
                Ez = (V(:,:,k-1)-V(:,:,k) + V(:,:,k)-V(:,:,k))/2;
                E(:,:,k) = sqrt(Exyx.^2+Exyy.^2+Ez.^2);
            end
            % max or mean?
            Field = max(nonzeros(elements_within_radius_3D(E,[mpx,mpy,mpz],[core_radius, core_radius+1])));
            density = perm * Field;
            Q = density * 4 * pi * core_radius * core_radius;
            numericalCap = Q/Volt;
            analyticalCap = (4*pi*perm)/((1/core_radius)-(1/outer_radius));
            error(err_start+z/check_every) = analyticalCap - numericalCap;
        end
        %}
    end
    E = zeros(Nx,Ny,Nz);
    % Can be optimised by making it only calculate E within the desired
    % radius
    for k = 2:Nz-1
        [Exyx, Exyy] = gradient(V(:,:,k)');%-gradient(permute(conj(V),[1,2,3]));
        Exyx = -Exyx;
        Exyy = -Exyy;
        Ez = (V(:,:,k-1)-V(:,:,k+1))/2;
        E(:,:,k) = sqrt(Exyx.^2+Exyy.^2+Ez.^2);
    end
    % max or mean?
    % use the outer radius this time - more smooth
    %Field = max(nonzeros(elements_within_radius_3D(E,[mpx,mpy,mpz],[core_radius, core_radius+1])));
    Field = max(nonzeros(elements_within_radius_3D(E,[mpx,mpy,mpz],[outer_radius-1, outer_radius+1])));
    density = perm * Field;
    Q = density * 4 * pi * outer_radius * outer_radius;
    %((outer_radius/core_radius)/0.4)-2
    numericalCap(1,col) = Q/Volt;
    analyticalCap(1,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius));
    errors(1,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius)) - Q/Volt;
    
    Field = mean(nonzeros(elements_within_radius_3D(E,[mpx,mpy,mpz],[outer_radius-1, outer_radius])));
    density = perm * Field;
    Q = density * 4 * pi * outer_radius * outer_radius;
    numericalCap(2,col) = Q/Volt;
    analyticalCap(2,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius));
    errors(2,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius)) - Q/Volt;
    
    Field = max(nonzeros(elements_within_radius_3D(E,[mpx,mpy,mpz],[core_radius-1, core_radius+1])));
    density = (perm * Field);
    Q = density * 4 * pi * core_radius * core_radius;
    %((outer_radius/core_radius)/0.4)-2
    numericalCap(3,col) = Q/Volt;
    analyticalCap(3,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius));
    errors(3,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius)) - Q/Volt;
    
    Field = mean(nonzeros(elements_within_radius_3D(E,[mpx,mpy,mpz],[core_radius, core_radius+1])));
    density = perm * Field;
    Q = density * 4 * pi * core_radius * core_radius;
    numericalCap(4,col) = Q/Volt;
    analyticalCap(4,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius));
    errors(4,col) = (4*pi*perm)/((1/core_radius)-(1/outer_radius)) - Q/Volt;
end

%% 3D plotting
% we can only really plot a slice of the data, offset to see up or down
offset = -60;

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
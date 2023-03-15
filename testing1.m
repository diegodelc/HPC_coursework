clear
close all


%% Setting up the grid
dx = 1; %dx = dy
Nx = 50;
Ny = 30;

u_grid_nopadding = zeros(Ny,Nx);
v_grid_nopadding = zeros(Ny,Nx);
h_grid_nopadding = zeros(Ny,Nx);

%Apply initial condition
%u and v initialsed at zero
for x = 1:Nx
    h_grid_nopadding(:,x) = 10 + exp(-(x-50)^2/25);
end


%Padding with the boundary conditions
stencil_length = 7;
u_grid_padded = applyPeriodicBcs(u_grid_nopadding,stencil_length);
v_grid_padded = applyPeriodicBcs(v_grid_nopadding,stencil_length);
h_grid_padded = applyPeriodicBcs(h_grid_nopadding,stencil_length);

%Calculating ks

%solving system in time


%plotting
[X,Y] = meshgrid(1:Ny,1:Nx)
figure()
plot3(X,Y,h_grid_nopadding)

%plotting padded grid
[X_pad,Y_pad] = meshgrid(1:Ny+6,1:Nx+6)
figure()
plot3(X_pad,Y_pad,h_grid_padded)


%% Functions
function grid_padded = applyPeriodicBcs(grid_nopadding,stencil_length)

    pad = (stencil_length-1)/2;

    [Ny,Nx] = size(grid_nopadding);
    Nx_padded = Nx + pad * 2;
    Ny_padded = Ny + pad * 2;

    grid_padded = zeros(Ny_padded,Nx_padded);
    
    grid_padded(4:36-3,4:56-3) = grid_nopadding;

   
    %top_rows
    grid_padded(1:pad,pad+1:Nx_padded-pad) = grid_nopadding(Ny-pad+1:Ny,:);

    %bottom rows
    grid_padded(Ny_padded-pad+1:Ny_padded,pad+1:Nx_padded-pad) = grid_nopadding(1:pad,:);

    %left columns
    grid_padded(pad+1:Ny_padded-pad,1:pad) = grid_nopadding(:,Nx-pad+1:Nx);

    %right columns
    grid_padded(pad+1:Ny_padded-pad,Nx_padded-pad+1:Nx_padded) = grid_nopadding(:,1:pad);

end
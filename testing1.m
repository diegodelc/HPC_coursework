clear
close all


%% Setting up the grid
dx = 1; %dx = dy
global Ny 
global Nx
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
stencil = (1/dx) * [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];
stencil_length = length(stencil);
u_grid_padded = applyPeriodicBcs(u_grid_nopadding,stencil_length);
v_grid_padded = applyPeriodicBcs(v_grid_nopadding,stencil_length);
h_grid_padded = applyPeriodicBcs(h_grid_nopadding,stencil_length);

%Calculating derivatives
dhdx = deriveX(h_grid_padded,stencil,stencil_length);
dhdy = deriveY(h_grid_padded,stencil,stencil_length);

dudx = deriveX(u_grid_padded,stencil,stencil_length);
dudy = deriveY(u_grid_padded,stencil,stencil_length);

dvdx = deriveX(v_grid_padded,stencil,stencil_length);
dvdy = deriveY(v_grid_padded,stencil,stencil_length);

dhudx = deriveX(h_grid_padded.*u_grid_padded,stencil,stencil_length);

dhvdy = deriveY(h_grid_padded.*v_grid_padded,stencil,stencil_length);

%calculating the ks

%solving system in time 


%plotting
[X,Y] = meshgrid(1:Ny,1:Nx);
figure()
plot3(X,Y,h_grid_nopadding)

%plotting padded grid
[X_pad,Y_pad] = meshgrid(1:Ny+6,1:Nx+6);
figure()
plot3(X_pad,Y_pad,h_grid_padded)


%% Functions
function grid_padded = applyPeriodicBcs(grid_nopadding,stencil_length)

    pad = (stencil_length-1)/2;

    [Ny,Nx] = size(grid_nopadding);
    Nx_padded = Nx + pad * 2;
    Ny_padded = Ny + pad * 2;

    grid_padded = zeros(Ny_padded,Nx_padded);
    
    grid_padded(pad+1:Ny_padded-pad,pad+1:Nx_padded-pad) = grid_nopadding;

   
    %top_rows
    grid_padded(1:pad,pad+1:Nx_padded-pad) = grid_nopadding(Ny-pad+1:Ny,:);

    %bottom rows
    grid_padded(Ny_padded-pad+1:Ny_padded,pad+1:Nx_padded-pad) = grid_nopadding(1:pad,:);

    %left columns
    grid_padded(pad+1:Ny_padded-pad,1:pad) = grid_nopadding(:,Nx-pad+1:Nx);

    %right columns
    grid_padded(pad+1:Ny_padded-pad,Nx_padded-pad+1:Nx_padded) = grid_nopadding(:,1:pad);

end

function d_wrt_x = deriveX(padded_data,stencil,stencil_length)
    global Ny;
    global Nx;
    d_wrt_x  = zeros(Ny,Nx);
    for xpos = 1:Nx %iterating over the columns, x
        for row = 1:Ny
           d_wrt_x(row,xpos) = sum(stencil.*(padded_data(row,xpos:xpos+stencil_length-1)));
        end
    end

end

function d_wrt_y = deriveY(padded_data,stencil,stencil_length)
    global Ny;
    global Nx;
    d_wrt_y  = zeros(Ny,Nx);
    for ypos = 1:Ny %iterating over the rows, y
        for col = 1:Ny
           d_wrt_y(ypos,col) = sum(stencil'.*(padded_data(ypos:ypos+stencil_length-1,col)));
        end
    end

end
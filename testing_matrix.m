clear
close all


%% The stuff
tmax = 5;
% Define grid and apply initial conditions

dx = 1; %dx = dy
global Ny 
global Nx
Nx = 100;
Ny = 100;


% myZeros = zeros(Ny,Nx);
u = zeros(Ny,Nx);
v = zeros(Ny,Nx);
% h = zeros(Ny,Nx);
global g
g = 9.81;

%initial conditions functions
which_test_case = 1;

H = 10;
fun1 = @(x,y) H + exp(-(x-50).^2./25);
fun2 = @(x,y) H + exp(-(y-50).^2./25);
fun3 = @(x,y) H + exp(-((x-50).^2 + (y-50).^2)./25);
fun4 = @(x,y) H + exp(-((x-25).^2+(y-25).^2)./25) + exp(-((x-75).^2+(y-75).^2)./25);

switch which_test_case
    case 1
        initial_cond = fun1;
    case 2
        initial_cond = fun2;
    case 3
        initial_cond = fun3;
    case 4 
        initial_cond = fun4;
end
[X,Y] = meshgrid(1:Nx,1:Ny);

h = initial_cond(X,Y);

%% Time thing
global stencil
stencil = (1/dx) * [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];

t = 0;
figure()
surf(X,Y,h)
title("start")
hold off
u = reshape(u,[Nx*Ny,1]);
v = reshape(v,[Nx*Ny,1]);
h = reshape(h,[Nx*Ny,1]);
yn = [u;v;h];
dt = 0.1;
counter = 0;
while t <= tmax
    
    
    
    
    
    k1 = calcFBLAS(yn);
    k2 = calcFBLAS(yn + dt*k1/2);
    k3 = calcFBLAS(yn + dt*k2/2);
    k4 = calcFBLAS(yn + dt*k3);

    yn = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)*dt; %updating data next timestep
    
%     size(yn(2*Ny+1:3*Ny,:))
%     clf    
%     surf(X,Y,yn(2*Ny+1:3*Ny,:),'EdgeColor',"None")
%     axis equal
%     drawnow

%     pause(0.1)
    
%     t = t+dt;
%     counter = counter +1;
end
%% Functions
%calculate derivatives and stack them
function F = calcFBLAS(yn)
    
    global Ny;
    global Nx;
    u_vect = yn(1:Nx*Ny,1);
    v_vect = yn(Nx*Ny+1:2*Nx*Ny,1);
    h_vect = yn(2*Nx*Ny+1:end,1);
    
%     size(u_vect)
%     size(v_vect)
%     size(h_vect)
    
    u = reshape(u_vect,[Nx,Ny]);
    v = reshape(v_vect,[Nx,Ny]);
    h = reshape(h_vect,[Nx,Ny]);
    %x
    dudx = derX(u);
    dvdx = derX(v);
    dhdx = derX(h);

    dudx = reshape(dudx,[Nx*Ny,1]);
    dvdx = reshape(dvdx,[Nx*Ny,1]);
    dhdx = reshape(dhdx,[Nx*Ny,1]);

    %y
    dudy = derY(u);
    dvdy = derY(v);
    dhdy = derY(h);

    dudy = reshape(dudy,[Nx*Ny,1]);
    dvdy = reshape(dvdy,[Nx*Ny,1]);
    dhdy = reshape(dhdy,[Nx*Ny,1]);


    ddx = zeros(3*Nx*Ny,1);
    ddy = zeros(3*Nx*Ny,1);

    counter = 1;
    for i = 1:Nx*Ny
        ddx(counter) = dudx(i);
        ddx(counter+1) = dvdx(i);
        ddx(counter+2) = dhdx(i);

        ddy(counter) = dudy(i);
        ddy(counter+1) = dvdy(i);
        ddy(counter+2) = dhdy(i);

        counter = counter + 3;
    end
    
    %reshaping u,v,h
    

    %populating matrices          < F = A * ddx + B*ddy >
    % a
    A = zeros(3*Nx*Ny,3*Nx*Ny);
    for i=1:Nx*Ny
        diag = (i-1)*3 + 1;
        A(diag,diag) = -u_vect(i);
        A(diag+1,diag+1) = -u_vect(i);
        A(diag+2,diag+2) = -u_vect(i);

        A(diag,diag+2) = -9.81;
        A(diag+2,diag) = h_vect(i);    
    end


    % b
    B = zeros(3*Nx*Ny,3*Nx*Ny);
    for i=1:Nx*Ny
        diag = (i-1)*3 + 1;
        B(diag,diag) = -v_vect(i);
        B(diag+1,diag+1) = -v_vect(i);
        B(diag+2,diag+2) = -v_vect(i);

        B(diag+1,diag+2) = -9.81;
        B(diag+2,diag+1) = -h_vect(i);
    end
    
    F = A*ddx + B*ddy;
end

%derivative in y
function derivative_grid = derY(data)
    global Ny;
    global Nx;
    global stencil;

    derivative_grid = zeros(Ny,Nx);
    for col = 1:Nx
        for row = 1:Ny
           if 4<=row && row<=Ny-3
           derivative_grid(row,col) = sum(stencil.*...
                   [ data(row-3,col) , data(row-2,col), data(row-1,col),...
                     data(row,col) , data(row+1,col), data(row+2,col), ...
                     data(row+3,col)]);
           elseif row == 1
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(Ny-2,col) , data(Ny-1,col), data(Ny,col),...
                     data(row,col) , data(row+1,col), data(row+2,col), ...
                     data(row+3,col)]);
           elseif row == 2
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(Ny-1,col) , data(Ny,col), data(row-1,col),...
                     data(row,col) , data(row+1,col), data(row+2,col), ...
                     data(row+3,col)]);
           elseif row == 3
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(Ny,col) , data(row-2,col), data(row-1,col),...
                     data(row,col) , data(row+1,col), data(row+2,col), ...
                     data(row+3,col)]);
           elseif row == Ny
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(row-3,col) , data(row-2,col), data(row-1,col),...
                     data(row,col) , data(1,col), data(2,col), ...
                     data(3,col)]);
           elseif row == Ny-1
              derivative_grid(row,col) = sum(stencil.*...
                   [ data(row-3,col) , data(row-2,col), data(row-1,col),...
                     data(row,col) , data(row+1,col), data(1,col), ...
                     data(2,col)]);
           elseif row == Ny-2
            derivative_grid(row,col) = sum(stencil.*...
                   [ data(row-3,col) , data(row-2,col), data(row-1,col),...
                     data(row,col) , data(row+1,col), data(row+2,col), ...
                     data(1,col)]);
           else
               disp("Error calculating derivative wrt y")
           end
        end
    end
end

%derivative in x
function derivative_grid = derX(data)
    global Ny;
    global Nx;
    global stencil;

    derivative_grid = zeros(Ny,Nx);
    for col = 1:Nx
        for row = 1:Ny
           if 4<=col && col<=Nx-3
           derivative_grid(row,col) = sum(stencil.*...
                   [ data(row,col-3) , data(row,col-2), data(row,col-1),...
                     data(row,col) , data(row,col+1), data(row,col+2), ...
                     data(row,col+3)]);
           elseif col == 1
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(row,Nx-2) , data(row,Nx-1), data(row,Nx),...
                     data(row,col) , data(row,col+1), data(row,col+2), ...
                     data(row,col+3)]);
           elseif col == 2
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(row,Nx-1) , data(row,Nx), data(row,col-1),...
                     data(row,col) , data(row,col+1), data(row,col+2), ...
                     data(row,col+3)]);
           elseif col == 3
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(row,Nx) , data(row,col-2), data(row,col-1),...
                     data(row,col) , data(row,col+1), data(row,col+2), ...
                     data(row,col+3)]);
           elseif col == Nx
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(row,col-3) , data(row,col-2), data(row,col-1),...
                     data(row,col) , data(row,1), data(row,2), ...
                     data(row,3)]);
           elseif col == Nx-1
               derivative_grid(row,col) = sum(stencil.*...
                   [ data(row,col-3) , data(row,col-2), data(row,col-1),...
                     data(row,col) , data(row,col+1), data(row,1), ...
                     data(row,2)]);
           elseif col == Nx-2
            derivative_grid(row,col) = sum(stencil.*...
                   [ data(row,col-3) , data(row,col-2), data(row,col-1),...
                     data(row,col) , data(row,col+1), data(row,col+2), ...
                     data(row,1)]);
           else
               disp("Error calculating derivative wrt x")
           end
        end
    end
end




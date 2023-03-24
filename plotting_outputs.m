clc
clear



% read data
path = "./cw_workspace/serial_implementation/output.txt";
data = importdata(path);

X = data(:,1);
Y = data(:,2);

u = data(:,3);
v = data(:,4);
h = data(:,5);

Nx = max(X)+1;
Ny = max(Y)+1;



H = zeros(Nx,Ny);
for x = 0:Nx-1
    
    for y = 0:Ny-1
        ind = x*Ny+y;
        H(x+1,y+1) = data(ind+1,5);
    end
end

[X,Y] = meshgrid(1:Nx,1:Ny); 

X = X';
Y = Y';



figure()
contourf(X,Y,H);
xlabel('x')
ylabel('y')
colorbar;

%figure()
%surf(X,Y,H);



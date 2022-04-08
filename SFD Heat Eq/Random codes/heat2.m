L    = 1; % Distance
N    = 100; % Number of grid points
nu   = 0.2; 
mesh = linspace(0,L,N);
dx   = mesh(2) - mesh(1); 
for i=1:N-2
    x(i)=(i-1)*dx; %since Matlab begins at 1, x_0 must be i-1
end
tsteps = 1000; % Number of time steps
tend = 2; 
t = linspace(0, tend, tsteps);
dt=1/500;
s = nu*dt/dx^2;
for i = 1:N-2
    for j = 1:N-2
        if j <= i - 2
            matrix(i,j) = 0;
        elseif j >= i + 2
            matrix(i,j) = 0;
        else
            if mod(i,2) ~= 0
                if mod(j,2) ~= 0
                    matrix(i,j) = 1+2*s;
                else
                    matrix(i,j) = -s;
                end
            else
                if mod(j,2) ~= 0
                    matrix(i,j) = -s;
                else
                    matrix(i,j) = 1+2*s;
                end
            end
        end
    end
matrix

plot(x, matrix(:,98))
end
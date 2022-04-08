N=10;
x = linspace(0, 1, N+1);
dx = x(2) - x(1);
global A B D

A = eye(N-1);
A = A *(-2 / dx^2);
A = A + 1/dx^2 * diag(ones(N-2,1), 1);
A = A + 1/dx^2 * diag(ones(N-2,1), -1);

B=zeros(N-1,1);

D=1/4;

u = exp(-(x-0.5).^2 / .25^2);
plot(x, u);
 
u0=u(2:end-1);
u0=u0';


[t, u] = ode45(@ddt_heat, [0, 1], u0);

plot(u)

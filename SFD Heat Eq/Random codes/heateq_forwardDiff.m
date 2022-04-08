% Initialize Data
% Length of Rod, Time Interval
% Number of Points in Space, Number of Time Steps
L=1;
T=0.1;
k=1;
N=10;
M=50;
dx=L/N;
dt=T/M;
alpha=k*dt/dx^2;

% Position
for i=1:N+1
x(i)=(i-1)*dx;
end
% Initial Condition
for i=1:N+1
u0(i)=sin(pi*x(i));
end
% Partial Difference Equation (Numerical Scheme)
for j=1:M
    for i=2:N
        u1(i)=u0(i)+alpha*(u0(i+1)-2*u0(i)+u0(i-1));
    end
    u1(1)=0;
    u1(N+1)=0;
    u0=u1;
end
% Plot solution
plot(x, u1);
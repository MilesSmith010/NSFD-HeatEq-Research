% Heat equation  in 1D 
% The PDE for 1D heat equation is Ut=Uxx, 0=<t,0=<x=<L
% Initial condions are U(0,t)=a(t);U(L,t)=b(t)
% the boundary condition is U(x,0)=g(x)
%  u(t,x) is the solution matrix.
% the finite linear heat equation is  solved is....
% -u(i-1,j)=alpha*u(i,j-1)-[1+2*alpha]*u(i,j)+alpha*u(i,j+1)...(1)
%alpha=dx/dt^2. dx,dt are  finite division for x and t.
% t is columnwise
%x is rowwise dealt in this code
%suggestions and discussions are welcome.
% All the best salorajan@gmail.com
clear;
clc;
clf;
% Input Data
L=1;tf=.01;
a = @(t) 0;         %2*t;       %boundary condition x=0
b = @(t) 0;           %2*t+L;     %boundary condition x=L
g = @(x) sin(pi.*x);  %x.^2;      %initial condition, t=0
dt=0.01;
dx=0.1;
alpha=dt/dx^2;
%The analytical solution is
ua=@(x,t) x.^2+2*t;
%------------------------
nx=L/dx;
nt=tf/dt;
t1=0:dt:tf;
x1=0:dx:L;
u=zeros(nt+1,nx+1);
u(:,1)=a(t1)';
u(:,end)=b(t1');
u(1,:)=g(x1);
%------------------------
%We observe the equation (1) is in form_ A*u(t,x)=-u(t-1,x)
%A is (nx+1,nx+1) size square natrix
A=eye(nx+1);
for i=2:nx
    A(i,i-1:i+1)=[alpha -1-2*alpha alpha];
end
% u(t,x) is deternuned sequentially
for i=2:nt+1
    B=-u(i-1,:)';
        B(1)=u(i,1);
        B(end)=u(i,end);
        u1 =(A\B)'; 
    u(i,2:nx) =u1(2:nx); 
end
%---------------------
%Plotting the results.
figure(1)
X1=kron(ones(nt+1,1),x1);
T1=kron(ones(1,nx+1),t1');
surf(X1,T1,u)
view(30,40);
title('Heat Equation 3D solution view')
xlabel('x')
ylabel('t')
zlabel('u')
figure(2)
for i=1:nt+1
plot(x1,u(i,:),'linewidth',2);
hold on;
end
title('Heat Equation 2D views for different times')
xlabel('x')
%----------------
%Juat change input data to solve 1D Heat equation. Thank you.
ylabel('u')
legend(num2str(t1'),'Location','northwest');

L=1;
tf=1;


a=@(t) 0*t;%2*t; %boundary condition x=0
b=@(t) 0;%2*t+L;%boundary condition x=L
g=@(x) sin(pi.*x) %x.^2; %initial condition, t=0


dt=0.01;
dx=0.1;


alpha=dt/dx^2; 

nx=L/dx
nt=tf/dt;

t1=0:dt:tf;
x1=0:dx:L;

u=zeros(nt+1,nx+1);

u(:,1)=a(t1')
u(:,end)=b(t1')
u(1,:)=g(x1)

%We observe the equation (1) is in form_ A*u(t,x)=-u(t-1,x)
%A is (nx+1,nx+1) size square natrix
A=eye(nx+1)
for i=2:nx
    A(i,i-1:i+1)=[alpha -1-2*alpha alpha];
end

A
%{
% u(t,x) is deternuned sequentially
for i=2:nt+1
    B=-u(i-1,:)';
        B(1)=u(i,1);
        B(end)=u(i,end);
        u1 =(A\B)'; 
    u(i,2:nx) =u1(2:nx); 
end
u
B

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
%}
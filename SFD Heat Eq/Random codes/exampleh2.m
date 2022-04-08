L = 1;
T = 1;
t_0 = 0;
x_0 = 0;
M=2^3;
K=100;

dt = T/K;
dx = L/M;

a = @(t) 0;         %2*t;       %boundary condition x=0
b = @(t) 0;           %2*t+L;     %boundary condition x=L
g = @(x) sin(pi.*x)  %x.^2      %initial condition, t=0

alpha=dt/dx^2; 

tk=0:dt:T;
xm=0:dx:L


%u (time = column, space = row)
u=zeros(K+1,M+1);

u(:,1)=a(tk'); %input BCs
u(:,end)=b(tk');
u(1,:)=g(xm) %input ICs

%We observe the equation (1) is in form_ A*u(t,x)=-u(t-1,x)
%A is (nx+1,nx+1) size square natrix
A=eye(M+1);
for i=2:M
    A(i,i-1:i+1)=[alpha -1-2*alpha alpha];
end

% u(t,x) is deternuned sequentially
for j=2:K+1
    B=-u(j-1,:)';
        B(1)=u(j,1);
        B(end)=u(j,end);
        u1 =(A\B)'; 
    u(j,2:M) =u1(2:M); 
end

%{
figure(1)
X1=kron(ones(K+1,1),xm);
T1=kron(ones(1,M+1),tk');
surf(X1,T1,u)
view(30,40);
title('Heat Equation 3D solution view')
xlabel('x')
ylabel('t')
zlabel('u')

%}
figure(2)
for i=1:M+1
    plot(xm,u(i,:),'linewidth',2);
    hold on;
end
title('Heat Equation 2D views for different times')
xlabel('x')


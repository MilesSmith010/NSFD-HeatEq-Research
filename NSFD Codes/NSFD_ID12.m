%CASE II ID#12

% Initialize
clc
clear
format long
po = 3;

L=1; %x length
T=1; %time to elapse
D=1/2; %Diffusivity coeff
M=2^po; %various different space step iterations

x_0=0;
t_0=0;

K=2^(2*po+1);
dx=L/M;
dt=T/K;

for j=1:K+1
    t(j)=(j-1)*dt;
end %time steps for mesh
t
t = 0:dt:T


R=D*dt/dx^2; % 
if R >= 1/2
    disp('R must be less than 1/2. Try lowering value of M')
    return
end


% Numerical Soln Loop

% Position Mesh Points
for i=1:M+1
    x(i)=(i-1)*dx; %since Matlab begins at 1, x_0 must be i-1
end

x = 0:dx:L; %for -L to L, change M to 2*M

a = @(t) 0;         %2*t;       %boundary condition x=0
b = @(t) 0;           %2*t+L;     %boundary condition x=L
g = @(x) sin(pi.*dx);  %x.^2;      %initial condition, t=0


u=zeros(K+1, M+1);

u(:,1)=a(t'); %input BCs
u(:,end)=b(t');
u(1,:)=g(x); %input ICs



% Forward Difference Scheme
for j=1:K
    for i=2:M
        u(j+1, i) = (1-2*R)*u(j,i) + (R*(u(j,i+1) - u(j,i-1))^2)/(2*u(j,i));
    end
end
  
    

%loop for exact solution
u_true=zeros(K+1, M+1);

for j=1:K+1
    for i=1:M+1
        u_true(j,i) = sin(pi.*x(i))*exp(-D*pi^2*t(j));       %note, the exp coul 
        
        %u(j+1,i)=(1-2R)*u(j,i) + R*(u(j,i+1) + u(j, i-1)) as alt equation
    end
end

%{
figure(1)
for j=1:K+1
    plot(x,u_true(j,:),'linewidth',2);
    hold on;
end
title('True Heat Equation 2D views for different times')
xlabel('x')
ylabel('u_truee')

figure(2)
for j=1:K+1
    plot(x,u(j,:),'linewidth',2);
    hold on;
end
title('Numerical Heat Equation 2D views for different times')
xlabel('x')
ylabel('u')
%}
figure(3)
X1=kron(ones(K+1,1),x);
T1=kron(ones(1,M+1),t');
surf(X1,T1,u)
view(30,40);
title('Numerical Heat Equation 3D solution view')
xlabel('x')
ylabel('t')
zlabel('u')

figure(4)
X1=kron(ones(K+1,1),x);
T1=kron(ones(1,M+1),t');
surf(X1,T1,u_true)
view(30,40);
title('True Heat Equation 3D solution view')
xlabel('x')
ylabel('t')
zlabel('u')
%{
figure(5)
plot(x,u(K,:),'linewidth',2);
hold on;
plot(x,u_true(K,:),'linewidth',2);
%}

% Errors and norms
    
u1 = u(end-1,:)
u1_true = u_true(end-1,:)


    L_inf = max(max(abs(u1-u1_true))); %NEED COMPUTE VECTOR NORM NOT MATRIX NORM, SO USE LAST ROW, COMPARE AT T = 1

    L2 = norm(u1-u1_true, 1); %FIX THIS

    L1  = norm(u1-u1_true, 2);

    L1 = L1';
    L2 = L2';
    L_inf = L_inf';
    
    x(2) - x(1)
    results = table(K, M, L1,L2,L_inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  NSFD Schemes ID #21 using Forward Difference Method for Heat Equation
  Explicit Equation:   
         U(j+1,k) = (1-2R)*U(j,k) + R(U(j,k+1) + U(j,k-1))

  Output: u := numerical approximation of u, the heat equation
          u_true := the true solution to the heat equation

   (CASE IV)
   ID#20 Unity Approximation:
         1 = (U(j,i) + U(j,i+1) + U(j,i-1))/(2U(j,i+1)

         U(j+1,k) = ((1-2*R)*U(j,i)*K)/(3U(j,i-1)) + ((1+R)(U(j,i-1)
         +K))/(3), K = u(j,i+1) + u(j,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables and mesh space
clc
clear
format long

L=1;                %x length, could be thought of as rod length
T=1;                %total t, typically 1
D=1/2;              %Diffusivity coefficient in (0,1)

po = 2;             %controls time and space steps
M=2^po;             %various different space step iterations
K=2^(2*po+1);       %this is set up such that K and M meet the CFL conditions for stability

x_0=0;              %initial x and t values
t_0=0;

dx=L/M;             %dx and dt are the evenly-spaced mesh points
dt=T/K;

R=D*dt/dx^2;        %CFL says R<1/2
if R >= 1/2
    disp('R must be less than 1/2. Try lowering value of M')
    return
end



%%Creating the Mesh

{
for j=1:K+1
    t(j)=(j-1)*dt;          %There are multiple ways to create the time mesh, this is not the most efficient but is clear in setup
end
}
t = 0:dt:T;                  %This is the method I use

Position Mesh Points
{
for i=1:M+1
    x(i)=(i-1)*dx; %since Matlab begins at 1, x_0 must be i-1
    end
}
x = 0:dx:L; %for -L to L, change M to 2*M








%% Numerical Soln Loop


a = @(t) 0;         %2*t;       %boundary condition x=0
b = @(t) 0;           %2*t+L;     %boundary condition x=L
g = @(x) sin(pi.*x);  %x.^2;      %initial condition, t=0


u=zeros(K+1, M+1);

u(:,1)=a(t'); %input BCs
u(:,end)=b(t');
u(1,:)=g(x); %input ICs



Forward Difference Scheme
for j=1:K
    for i=2:M
        W = u(j,i+1) + u(j,i-1);
        u(j+1,i) = ((1-2*R)*(u(j,i)+W)+R*W)/(3) + (R*W^2)/(3*u(j,i));
    end
end
  
    

Computing the Exact Solution, u_true
u_true=zeros(K+1, M+1);

for j=1:K+1
    for i=1:M+1
        u_true(j,i) = sin(pi.*x(i))*exp(-D*pi^2*t(j));   
    end
end

%% Plots and Graphs

Figures 1 - 2 are 1-D representations of u. Since these aren't often
useful, they are commented out.
{
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
}

3-D Representations of u and u_true
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
{
figure(5)
plot(x,u(K,:),'linewidth',2);
hold on;
plot(x,u_true(K,:),'linewidth',2);
}

%% Errors and norms

We are taking the LAST row of u and u_true to compare them. This will
capture any truncation errors that have become exacerbated over iterations
    
u1 = u(end-1,:)
u1_true = u_true(end-1,:)


L_inf = max(max(abs(u1-u1_true))); 
L2 = norm(u1-u1_true, 1); 

L1  = norm(u1-u1_true, 2);

L1 = L1';
L2 = L2';
L_inf = L_inf';

x(2) - x(1)
results = table(K, M, L1,L2,L_inf)

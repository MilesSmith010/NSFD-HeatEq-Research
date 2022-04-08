%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NSFD Test using Backward Difference Method for Heat Equation
%   Explicit Equation:   
%          U(j+1,k) = A-1(u(j,i)+b)
%
%   Output: u := numerical approximation of u, the heat equation
%           u_true := the true solution to the heat equation
%
%    A := tridiagonal matrix of solutions for the time mesh
%    b = boundary values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialize variables and mesh space
clc
clear
format long

L=1;                %x length, could be thought of as rod length
T=1;                %total t, typically 1
D=1/2;              %Diffusivity coefficient in (0,1)

po = 4;             %controls time and space steps
M=2^po;             %various different space step iterations
K=2^po;       %this is set up such that K and M meet the CFL conditions for stability

x_0=0;              %initial x and t values
t_0=0;

dx=L/M;             %dx and dt are the evenly-spaced mesh points
dt=T/K;

t = 0:dt:T;         %This is the method I use
x = 0:dx:L;         %for -L to L, change M to 2*M


R=D*dt/dx^2;        %CFL condition unnecessary due to BWFD method unconditional convergnce

a = @(t) 0;         %2*t;       %boundary condition x=0
b = @(t) 0;           %2*t+L;     %boundary condition x=L
g = @(x) sin(pi.*x);  %x.^2;      %initial condition, t=0

u=zeros(K+1, M+1);

u(:,1)=a(t'); %input BCs
u(:,end)=b(t');
u(1,:)=g(x); %input ICs


d0 = 1+2*R;
d1 = -R;
%A = zeros(M-1,M-1)
A = diag(d0*ones(1,M+1)) + diag(d1*ones(1,M),1) + diag(d1*ones(1,M),-1)
A = inv(A);

b=zeros(1,M+1);
for j=1:K
    b(1) = R*u(j,1);
    b(M+1) = R*u(j,M+1);
    u(j+1,:) = (A*(u(j,:) - b)');
end

%Computing the Exact Solution, u_true
u_true=zeros(K+1, M+1);

for j=1:K+1
    for i=1:M+1
        u_true(j,i) = sin(pi.*x(i))*exp(-D*pi^2*t(j));   
    end
end


%3-D Representations of u and u_true
figure(3)
X1=kron(ones(K+1,1),x);
T1=kron(ones(1,M+1),t');
surf(X1,T1,u)
view(30,40);
title('Numerical Heat Equation 3D solution view')
xlabel('x')
ylabel('t')
zlabel('u')

    
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


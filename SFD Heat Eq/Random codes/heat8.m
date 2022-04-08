M = 10;
K = 10;
x = linspace(0, 1, M+1);
dx = x(2) - x(1);
t = linspace(0, 1, K+1);
dt = t(2) - t(1);
global A B D

A = eye(M-1);
A = A *(-2 / dx^2);
A = A + 1/dx^2 * diag(ones(M-2,1), 1);
A = A + 1/dx^2 * diag(ones(M-2,1), -1);

B=zeros(M-1,1);

D=1/4;



u_0 = exp(-(x-0.5).^2 / .25^2);
plot(x, u);
 
u0=u(2:end-1);
u(1,:)=u0

for j=1:K
    u(1,j+1)=
    for i=2:M
        u(i, j+1)=u(i,j)+R*(u(i+1,j)-2*u(i,j)+u(i-1,j));
    end
end


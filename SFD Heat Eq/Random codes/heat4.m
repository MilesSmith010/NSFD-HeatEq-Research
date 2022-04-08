%implicit euler for bernoulli 2
function [max_err,err_2, err_1,x] = bernoulli2_impliciteuler(h)
lambda =1;
% start computation
k=1;
k_bar=k;
x(1)=1/4; %matlab uses 1-indexing
u(1)=1/2;
while k_bar*h<=3
    k=k+1;
    k_bar=k;
    x=cat(2,x(1),zeros(1, k_bar-1));
    
    for i=2:k_bar
       a=lambda*h + 1;
       b=-lambda*h;
       c=-(u(i-1))^2;
       u(i)= (-b + sqrt(b^2 - 4*a*c))/(2*a);
       x(i)=(u(i))^2;
        % locates the index at which we lose positivity
    end  
end



x_true=zeros(1, k_bar);

% time interval
t(1)=0;
t = (t(1):h:h*(k_bar-1));    

% compute the true solution
for i = 1:k_bar
    x_true(i) = x_exact(t(i));
end


%{
 Calculate the MAXIMUM ERROR max_err between the exact solution x_true and 
 the numerical solution x
%}
max_err = max(abs(x-x_true));

%{
Calculate the L2 NORM ERROR err_2 between the exact solution x_true and the 
numerical solution x
%}
err_2 = sqrt(dot(x -x_true, x -x_true)*h);
err_2 = norm(x-x_true);

%{
Calculate the L1 NORM ERROR err_1 between the exact solution x_true and the 
numerical solution x
%}
err_1 = norm(x-x_true, 1);


end

%Script for calculating the residual.
function [ode_error]=ode_error_fun(parameters,x,tSpan)
%I assume opt is options and non are added therefore you supply matlab with
%an empty vector.
opt=[];
%I don't comprehend why r and b are inputs (as parameters of course) in the ODE
%script. Is it because the system is '2D' and therefore we can chose to fix
%two of the parameters at 1, or is it some similar assumption we make so that
%we can fix them?
r=1;
b=1;
%Adding the initial values
x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
%Turning them into a vector
X0 = [x1(1) x2(1) x3(1)];
%Calling the solver.
[~,yx] = ode45(@LSQodes,tSpan,X0,opt,parameters);
ode_error  = yx(:,:)-x(:,:);%Original script
%ode_error = Residual.*Residual;
%My version uses Chi-square as optimization
%ode_error = (x(:,:)-yx(:,:)).^2/var(x(:,:));
%Calculating the residual
end
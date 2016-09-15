function [Parameters]=LSQ2(P)
%LSQ takes one vector of 5 numbers as input.
%%These 7 numbers will be used as starting values for the parameters.
%Load the data
Data = dlmread('newdataset2withnoise.txt','\t');
%Lower bound that the algorithm is allowed to search
lowerBounds = [0.0001,0.0001,0.0001,0.0001,0.0001];
%Upper bound that the algorithm is allowed to search
uperBounds = [300,300,300,300,300];
%Values that LSQ start looking at. I'm not sure if this is redundant coding
%or not, could I just name the input 'StartParameters' and skip the 'p'
%annotation? 
%Adding options to the algorithm
options=optimset('Algorithm','levenberg-marquardt', 'MaxFunEvals',3000, 'MaxIter',2000,'TolFun',1e-5,'TolX', 1e-20,'Display', 'iter');
%Description of options:
%%MaxFunEvals = Maximum number of function evaluations allowed, a positive integer. The default is 100*numberOfVariables
%%MaxIter = Maximum number of iterations allowed, a positive integer. The default is 400
%%TolFun = function tolerance = ???
%%TolX = ???
%%Display = level of display
%%Iter = show output after each computation
%Adding the time data from the 'time column' in the data
tSpan = Data(:,1);
%Initial conditions extracted from the data
%Storing the parameters in L
Parameters = lsqnonlin(@(p) ode_error_fun(p,Data(:,[2 3 4]),tSpan),P,[],[],options);
plotLSQ2([Parameters 1 1]);
end
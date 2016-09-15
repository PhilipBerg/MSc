function [Pen,Pars,OptiPen,OptiPars]=LSQ(x,P,Lower,Upper,N)
%LSQ takes the amount of variables(x) and paramters(P) Lower- and Upper paramter value,
%and the N amount of evaluations. It returns the best found sum of residuals, and
%its corresponding paramters, and the paramters further optimized by lsqnonlin.
%Loading the data
Data = dlmread('dataset1withnoise.txt','\t');
%Adding options to the algorithm with levenberg-marquardt
options=optimset('Algorithm','levenberg-marquardt', 'MaxFunEvals',3000, 'MaxIter',2000,'TolFun',1e-20,'TolX', 1e-20,'Display', 'iter');
%Without levenberg-marquardt
%options=optimset('MaxFunEvals',3000, 'MaxIter',2000,'TolFun',1e-5,'TolX', 1e-20,'Display', 'iter');
%Description of options:
%%MaxFunEvals = Maximum number of function evaluations allowed, a positive integer. The default is 100*numberOfVariables
%%MaxIter = Maximum number of iterations allowed, a positive integer. The default is 400
%%TolFun = function tolerance = Termination tolerance on the function value
%%TolX = Termination tolerance on x.
%%Display = level of display
%%Iter = show output after each computation
%Adding the time data from the 'time column' in the data
tSpan = Data(:,1);
%Drawing N random random sets of p parameters between Upper- and Lower
%values, parameters are randomized unbiased.
ParameterSets = Lower+(Upper-Lower).*rand(N,P);
ii=1;
Penalties=zeros(N,x);
%Creating a loop that stores the total residual and parameters for each run
while (ii<=N),
    [res] = ode_error_fun(ParameterSets(ii,:),Data(:,[2 3 4]),tSpan);
    Penalties(ii,:)= sum(res);
    ii=ii+1;
end
%Concatenate Penalties and Parameters
[PenPar] = [Penalties ParameterSets];
%Looping step 3-7 10 times more then total amount of paramter sets
for BigLoop = 1:N*10
      %Create a (N/10)-by-(N/10) matrix of integers.
      [drawRow] = randi(N,round(N/10));
      %Take a tenth of the paramters randomly
      [RandomParameters] = PenPar(drawRow(:,1),x+1:end);
      %Pulling out one of the sets of parameters
      Pdraw= randi(length(RandomParameters),1);
      [Pn] = RandomParameters(Pdraw,:);
      %Removing it from the original set
      RandomParameters(Pdraw,:)=[];
      %Creating the new set of paramters
      Pnew=2.*mean(RandomParameters)-Pn;
      %%%%%%%%Remove negative or 0 values
      %If any value in Pnew<0 sum(debugger) has to be >0
      Debugger = (Pnew<0);
    while sum(Debugger) > 0,
        [RandomParameters] = PenPar(drawRow(:,1),x+1:end);
        Pdraw= randi(length(RandomParameters),1);
        Pn = RandomParameters(Pdraw,:);
        RandomParameters(Pdraw,:)=[];
        Pnew=2.*mean(RandomParameters)-Pn;
        Debugger = (Pnew<0);
    end
    %Calculating penalty for the new paramter vector
    [res2] = ode_error_fun(Pnew,Data(:,[2 3 4]),tSpan);
    PnewPenalty=sum(res2);
    SumResiduals=zeros(length(PenPar),1);
    %Lump all residuals together
    for loop = 1:length(PenPar)
        SumResiduals(loop,1) = sum(PenPar(loop,1:x));
    end
    %Checking if the new paramters have lower residual then the biggest
    %ones
    if sum(PnewPenalty) < max(SumResiduals)
        NewAddition=[PnewPenalty Pnew];
        [~,RowIndex] = max(SumResiduals);
        PenPar(RowIndex,:) = NewAddition;
    end
end
%Grab the lowest residuals
[~,RowDex]=min(SumResiduals);
%Grab the lowest residuals
[Pen]=PenPar(RowDex,1:x);
%Grab the best parameters
[Pars]=PenPar(RowDex,x+1:end);
%Run them through LSQnonlin
[OptiPars,OptiPen] = lsqnonlin(@(p) ode_error_fun(p,Data(:,[2 3 4]),tSpan),PenPar(RowDex,x+1:end),[],[],options);
%Plot the original parameters
figure(1);
plotLSQ(Pars);
%Plot the optimized parameters
figure(2);
plotLSQ(OptiPars);
end
function phase(p)
%Add empty figure
figure;
hold on;
%Vectors for different starting values
n=[0 10 20 30 40 60 80 100];
E=[100 80 60 40 30 20 10 0];
%Do as many simulations as starting values
i=length(n);
%Loop and plot simulations
while i>0,
  [~,dX] = ode45(@(t,x) LSQodes(t,x,[p 1 1]),[0 200],[200-n(i)-E(i), n(i), E(i)]);
   plot(dX(:,2),dX(:,3),... %drawing x,y
       'LineWidth',2);
   plot(dX(1,2),dX(1,3),'bo',...% starting points
       'MarkerSize',10); 
   plot(dX(end,2),dX(end,3),'ks',...% ending points
       'MarkerSize',10); 
   i=i-1;
end;
 %Creating even spaced matrixes with cordinates 0 to 50
 x2 = linspace(0,100,20);
 y2 = linspace(0,100,20);
 [x3,y3] = meshgrid(x2,y2);
    %Include Z if you want a 3D plot
    %Z = y3.^2 - x3.^2;
    %Include if you want a 3D plot
    %[X,Xm,Xstar] = surfnorm(Z);
 %Creating matrixes to hold the value of the derivative at t=0.
Xm = zeros(size(x3));
Xstar = zeros(size(x3));
%Adding the derivative of Xm and Xstar at t=0 (Yprime(1,2)&(1,3)) to the grid.
for Q = 1:numel(x3)
    dx = LSQodes(0,[200-x3(Q)-y3(Q), x3(Q), y3(Q)],[p 1 1]);
    Xm(Q) = dx(2);
    Xstar(Q) = dx(3);
end 
%Scaling the size of the arrows in the phase field
for Z = 1:numel(x3)
 Vmod = sqrt(Xm(Z)^2 + Xstar(Z)^2);
 Xm(Z) = Xm(Z)/Vmod;
 Xstar(Z) = Xstar(Z)/Vmod;
end
scale = 0.40;
%Calling quiver to plot the arrows
 quiver(x3,y3,Xm,Xstar, scale);
 %Adding nullclines
 %The nullclines of the current system are not very informative.
 %setting parameters
% k1 = p(1);
% k2 = p(2);
% alpha = p(3);
% beta = p(4);
% K = p(5);
% B = 1;
% R = 1;
 %Scaling the axis.
 %[t,Null] = LSQodes(t,x,[p 1 1]),[0 200],[100, 0, 0,0,0]);
 %plot(t,Null(:,2),... %drawing nullclines
%       'LineWidth',3);
% plot(t,Null(:,3),... %drawing nullclines
%       'LineWidth',3);
 xlabel('Xm');
 ylabel('Xstar');
 title('Phase field with simulations at different initial values');
axis tight equal;
 grid;
 hold off;
end
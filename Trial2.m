function [t, dydt1] = HematopoeiticSystemTrial
tspan = [0,4];
x0 = [1/100 1/100 1/100 1/100];

P0=0.7; 
h0=0.0000235;
r0=0.01818;
r1=0.08712;
r2=8.0217; 
d3=0.1;
p1=0.4783;
p2=0.4987; 
[t,y] = ode45(@(t,y)myODE(t,y,P0,h0,r0,r1,r2,d3,p1,p2),tspan,x0);

%--------Figures--------%
figure (1)
hold on
plot(t,y(:, 1), 'g')
plot (t, y(:,2),'y')
plot (t, y(:,3),'b')
plot (t, y(:,4), 'r')
end

function dydt = myODE(t,y,P0,h0,r0,r1,r2,d3,p1,p2);
p0=P0/(1+h0*y(1));
tmp1 = r0*y(1)*(p0-1);
tmp2 = 2*r0*y(1)*(1-p0)+r1*y(2)*(2*p1-1);
tmp3 = 2*r1*y(2)*(1-p1)+r2*y(2)*(2*p2-1);
tmp4 = 2*r2*y(3)*(1-p2)-d3*y(4);
dydt = [tmp1; tmp2; tmp3; tmp4];
end
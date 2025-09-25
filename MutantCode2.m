function [t,dydt1] = HematopoieticSystemModel

time = 0:10000;
%------%

x_0 = 17000;
y_0 = 0;
x_1 = 81800;
y_1 = 1;
x_2 = 357000;
y_2 = 0;
x_3 = 28687000;
y_3 = 0;
% p0_0=100;
%----------%

p = [x_0 y_0 x_1 y_1 x_2 y_2 x_3 y_3];
%--------------------------------------------------%
[t,dydt1] = ode45(@myODE,time,p);


%---Matrices-----%
% [m1, m2, m3, m4, m5, m6] = equillibriatest(dydt1(:,1), dydt1(:,2), dydt1(:,3), dydt1(:,4), dydt1(:,5), dydt1(:,6));
% save('m1.mat','m1')
% save('m2.mat','m2')
% save('m3.mat','m3')
% save('m4.mat','m4')
% save('m5.mat','m5')
% save('m6.mat','m6')


%% -----Individual Plots------ %
% figure (1)
% hold on
% subplot (1,3,1);
% plot(t,dydt1(:, 1), 'r');
% title ('HSC');
% xlabel('time')
% ylabel('cells')
% hold off
% 
% hold on
% subplot (1,3,2);
% plot (t, dydt1(:,3),'g');
% title ('ST-HSCs');
% xlabel('time')
% ylabel('cells')
% hold off
% 
% hold on
% subplot (1, 3, 3);
% plot (t, dydt1(:,5),'b');
% title ('MPPs');
% xlabel('time')
% ylabel('cells')
% hold off
% 
% %------------------%
figure (2)
hold on
plot(t,dydt1(:, 1), 'r','LineWidth',1.0);
plot (t, dydt1(:,3),'g','LineWidth',1.0);
plot (t, dydt1(:,5),'b','LineWidth',1.0);
plot (t, dydt1(:,7),'y','LineWidth',1.5);
plot(t,dydt1(:, 2), '--r','LineWidth',1.0);
plot (t, dydt1(:,4),'--g','LineWidth',1.0);
plot (t, dydt1(:,6),'--b','LineWidth',1.0);
plot (t, dydt1(:,8),'--y','LineWidth',1.5);
legend ('HSC' , 'ST-HSC', 'MPP', 'CLP/CMP');
xlabel('time (days)');
ylabel('cells');
ylim([10^-2 10^12]);
xlim([10^0 10^4]);
title ('Hematopoietic System Growth');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
hold off

%------------------%
% 
% figure (3)
% plot(t,dydt1(:, 2), 'g');
% hold on
% plot (t, dydt1(:,4),'b');
% plot (t, dydt1(:,6),'m');
% legend ('HSC' , 'ST-HSC', 'MPP');
% xlabel('time (days)');
% ylabel('mutant cells');
% title ('Hematopoietic System: Mutant Growth');
% set(gca, 'YScale', 'log');
% hold off
% ylim auto
% % 
% figure(4)
% plot (t, dydt1(:,1),'g');
% hold on
% plot (t, dydt1(:,2), 'k');
% title ('HSCs');
% xlabel('time');
% ylabel('cells');
% legend ('HSC', 'mutant');
% set(gca, 'YScale', 'log');
% hold off
% ylim auto

% figure (5)
% plot (t, dydt1(:,3),'b')
% hold on
% plot (t, dydt1(:,4),'k')
% title ('ST-HSCs');
% xlabel('time');
% ylabel('cells');
% legend ('ST-HSC', 'mutant');
% set(gca, 'YScale', 'log');
% hold off
% ylim auto
% 
% figure (6)
% plot (t, dydt1(:,5), 'm')
% hold on
% plot (t, dydt1(:,6),'k')
% title ('MPPs');
% xlabel('time');
% ylabel('cells');
% legend ('MPPs', 'mutants');
% set(gca, 'YScale', 'log');
% hold off 
% ylim auto


% figure(7)
% hold on
% plot (t, dydt1(:,5), 'k')
% title ('Self Renewal for HSCs');
% xlabel('time')
% ylabel('cells')
% hold off

end

function dydt = myODE(t,y, r);

global P0 P1 P2 h0 r0 r1 r2 d3 u1 s x0 u0 u2 h2 h1 
d3=0.1;
s=0.05;
x0=17000;
r0=0.01818;
r1=0.8712;
r2=8.0217;
P0=0.5;
P1=0.4783;
P2=0.4987;
u0=0;
u1=0;
u2=0;
%%weak control
h0= 10^-5;
h1= 10^-5;
h2= 10^-5;

% y(1) = x_0;
% y(2) = y_0;
% y(3) = x_1;
% y(4) = y_1;
% y(5) = x_2;
% y(6) = y_2;

c0=0.5+((h0*x0)/2);
c1=P1+((r0*x0*h1*P1)/(r1*(1-2*P1)));
c2=P2+((2*r0*x0*h2*P2*(1-P1))/(r2*(1-2*P1)*(1-2*P2)));

p0=c0./(1+h0.*(y(1)+y(2)));
p1=c1./(1+h1.*(y(3)+y(4)));
p2=c2./(1+h2.*(y(5)+y(6)));

p_0=(1+s)*c0/(1+h0*(y(1)+y(2)));
p_1=(1+s)*c1/(1+h1*(y(3)+y(4)));
p_2=(1+s)*c2/(1+h2*(y(5)+y(6)));


tmp1 = r0.*y(1).*p0.*(1-u0)-r0.*y(1).*(1-p0);
tmp2 = r0.*y(1).*p0.*u0-r0.*y(2).*(2.*p_0-1); 
%%
tmp3 = r0.*y(1).*(1-p0).*(2-u0)+r1.*y(3).*p1.*(1-u1)-r1.*y(3).*(1-p1);
tmp4 = r0.*y(1).*(1-p0).*u0+2.*r0.*y(2).*(1-p_0)+r1.*y(3).*p1.*u1+r1.*y(4).*(2*p_1-1);
%%
tmp5 =r1.*y(3).*(1-p1).*(2-u1)+r2.*y(5).*p2.*(1-u2)-r2.*y(5).*(1-p2);
tmp6 =r1.*y(3).*(1-p1).*u1+2.*r1.*y(4).*(1-p_1)+r2.*y(5).*p2.*u2+r2.*y(6).*(2*p_2-1);
%%
tmp7 =r2*y(5)*(1-p2)*(2-u2)-d3*y(7);
tmp8 =r2*y(5)*(1-p2)*u2+2*r2*y(6)*(1-p_2)-d3*y(8);

dydt = [tmp1; tmp2; tmp3; tmp4; tmp5; tmp6; tmp7; tmp8];

end

% function [matrix1, matrix2, matrix3, matrix4, matrix5, matrix6] = equillibriatest (in1, in2, in3, in4, in5, in6)
% matrix1 = in1;
% matrix2 = in2;
% matrix3 = in3;
% matrix4 = in4;
% matrix5 = in5;
% matrix6 = in6;
% end
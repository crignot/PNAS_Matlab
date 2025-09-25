function [t,dydt1] = HematopoieticSystemModel

time = 0:6000;
%------%

x_0 = 1;
x_1 = 3;
x_2 = 9;
x_3 = 50;
y_1 = 3;
% p0_0=100;
%----------%

p = [x_0 x_1 x_2 x_3 y_1];
%--------------------------------------------------%
[t,dydt1] = ode45(@myODE,time,p);


%---Matrices-----%
[m1, m2, m3, m4, m5] = equillibriatest(dydt1(:,1), dydt1(:,2), dydt1(:,3), dydt1(:,4), dydt1(:,5))
save('m1.mat','m1')
save('m2.mat','m2')
save('m3.mat','m3')
save('m4.mat','m4')
save('m5.mat','m5')


%% -----Individual Plots------ %
figure (1)
hold on
subplot (2,2,1);
plot(t,dydt1(:, 1), 'r');
title ('HSC');
xlabel('time')
ylabel('cells')
hold off

hold on
subplot (2,2,2);
plot (t, dydt1(:,2),'g');
title ('ST-HSCs');
xlabel('time')
ylabel('cells')
hold off

hold on
subplot (2, 2, 3);
plot (t, dydt1(:,3),'b');
title ('MPPs');
xlabel('time')
ylabel('cells')
hold off

hold on
subplot (2, 2, 4);
plot (t, dydt1(:,4), 'm');
title ('CLPs and CMPs');
xlabel('time')
ylabel('cells')
hold off

%------------------%
figure (2)
hold on
plot(t,dydt1(:, 1), 'r');
plot (t, dydt1(:,2),'g');
plot (t, dydt1(:,3),'b');
plot (t, dydt1(:,4), 'm');
legend ('HSC' , 'ST-HSC', 'MPP', 'CLP and CMP');
xlabel('time (days)');
ylabel('cells');
title ('Hematopoietic System Growth');
set(gca, 'YScale', 'log');
hold off

%------------------%

figure (3)
hold on
plot(t,dydt1(:, 1), 'r');
yline(1.70213e+04)
title ('HSC');
xlabel('time');
ylabel('cells');
hold off

figure(4)
hold on
plot (t, dydt1(:,2),'g');
yline(8.18424e+04);
plot (t, dydt1(:,5), 'b');
title ('ST-HSCs');
xlabel('time');
ylabel('cells');
legend ('ST-HSC', 'Equillibrium Level', 'mutants');
set(gca, 'YScale', 'log');
hold off

figure (5)
hold on
plot (t, dydt1(:,3),'b')
yline(3.56703e+05);
title ('MPPs');
xlabel('time');
ylabel('cells');
legend ('MPP', 'Equillibrium Level');
hold off

figure (6)
hold on
plot (t, dydt1(:,4), 'm')
title ('CLPs and CMPs');
yline(2.8687e+07)
xlabel('time');
ylabel('cells');
legend ('CLP and CMP', 'Equillibrium Level');
hold off 


% figure(7)
% hold on
% plot (t, dydt1(:,5), 'k')
% title ('Self Renewal for HSCs');
% xlabel('time')
% ylabel('cells')
% hold off

end

function dydt = myODE(t,y, r);

global P0 h0 r0 r1 r2 d3 p1 p2 u1 s
P0=0.7; 
h0=0.0000235;
r0=0.01818;
r1=0.08712;
r2=8.0217; 
d3=0.1;
p1=0.4783;
p2=0.4987; 
u1= 0.00001;
s=0.046;
%%

% y(1) = x_0;
% y(2) = x_1;
% y(3) = x_2;
% y(4) = x_3;

p_1=p1.*(1+s);
p0=P0./(1+h0.*y(1));
tmp1 = r0.*y(1).*(2.*p0-1);

tmp2 = r0.*y(1)+r1.*y(2).*p1.*(1-u1)-r1.*y(2).*(1-p1);

tmp3 = r0.*y(1)+r1.*y(2).*p1.*(1-u1)-r1.*y(2).*(1-p1);

tmp4 =2.*r2.*y(3).*(1-p2) - d3.*y(4);

tmp5 = r1.*y(2).*p1.*u1 +r1.*y(5).*(2*p_1-1);

dydt = [tmp1; tmp2; tmp3; tmp4; tmp5];

end

function [matrix1, matrix2, matrix3, matrix4, matrix5] = equillibriatest (in1, in2, in3, in4, in5)
global P0 h0 r0 r1 r2 d3 p1 p2 u1 s
matrix1 = in1;
matrix2 = in2;
matrix3 = in3;
matrix4 = in4;
matrix5 = in5;
end
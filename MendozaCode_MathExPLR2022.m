% Describe your ODE Here:

% This is an example
% ---- x'= (r/k)x(k-x)-C(xy/(1+ax)) -------%
% -----y'= -dy + b(xy/(1+ax))--------------%
% ----------------------------------------%


%In our model, we have 4 variables, see paper : 
% ---- P_0'=  (1/tau_0)(1-P_0- P_g)  +   (1/tau_g0) X_L(L) P_g  - (1/tau_o) X_L*(L)  X_G(G)P_0    -------%
% -----P_g'=  ...  ---------%  
% -----G'= ...    ---------%     
% -----L'=...------------%
% For convenienvce ,we can use x,y,z,r as 4 variables.

function exitflag = samplecodeforMathExPLR2022
%
% --------------------------------------------------%
% -----Set the initial values---------------------%
% --------------------------------------------------%

%time = [start_time: period :end_time]
time1 = [0:1:4]; 
time = [0:0.1:4]; 

% time = [0:5:20];

% -----  Initial values ---------------------%

%p_s=  1 *10^5; 


 %p_0 + p_g = 100  with scale 850 ;

 %  p_0 = ? 




 
 
 
 
 p_0 = 90/850  ;
 p_g = 100/850 - p_0;


%p_0= 0.2;
%p_g= 0.01;



G   = 1.7 ;
L   = 0.15;



%

%%%%  Real _ data %%%%%%


P_total = [100, 190, 320 ,440 580];  


G_real = [1.7 1.45  0.95 0.45 0.2 ];

L_real = [ 0.15 0.2 0.68 1.08 1.45];



%%%%%%%%%%%%%%%%%%%%%%%



% ----- Paremeters ---------------------%
% You can add the parameters here and give them values.

global tau_o tau_go p_star tau_og   tau_g lambda L_star Kk_o K_L K_G Gamma n_star  m_star  G_star 



tau_o =  1/log(2);  % known
tau_g =  1/ log(2);  % known

p_star = 850;     % known

tau_og = 1/24;   % known 

tau_go = 1;      % known

lambda = 1;    % known

L_star = 10;      % known

Kk_o = 1.2;       %  known

K_L = 5;         %known

K_G = 2 ;      % known

Gamma = 100 ; % known

n_star = 1;% known

m_star = 1;% known

G_star = 0.5;% known





% Put all of the Initialization variables together, like this 

p = [p_0 p_g G L ];
%--------------------------------------------------%
%----using ODE45 to solve the equations  -----%
%--------------------------------------------------%
[t,dydt1] = ode23(@Mondoza,time,[p(1),p(2),p(3),p(4)],[],p);



%---- for 4 variables--- %%
%[t,dydt] = ode45(@yourfunction,time,[p(1),p(2)],[],p);
% The results dydt has include 4 variables , for x=>  dydt(:,1)   y=>% dydt(:,2) ans so on

%Plot figures 


Total =   (dydt1(:,1) + dydt1(:,2)) * 850 ; 
% subplot(3,1,1)
% hold on
% 
% plot(time, Total,'-o')
% plot (time ,P_total,'-*')
% % axis([0 4  0 700]);
% hold off
% 
% 
% subplot(3,1,2)
% 
% hold on
% plot(time, dydt(:,3),'-o')
% plot (time ,G_real,'-*')
% % axis([0 4  0 1.7]);
% hold off
% 
% 
% subplot(3,1,3)
% 
% hold on
% plot(time, dydt(:,4),'-o')
% plot (time ,L_real,'-*')
% % axis([0 4  0 1.5]);
% hold off



figure(1)
hold on

plot(time, Total,'-')
%plot(time, dydt1(:,1)*850)
%plot(time, dydt1(:,2)*850)

plot (time1 ,P_total,'v')
set(gca,'ytick',0:200:630)
set(gca,'xtick',0:1:4)
 axis([0 4  0 630]);
%legend ('Total cells' , '# of P_o cells ', '# of P_g cells')
xlabel('days in culture')
ylabel('cells/ml(x10^3)')


hold off


figure(2)

hold on
plot(time, dydt1(:,3),'-')
plot (time1 ,G_real,'v')
%legend ('simulation', 'Data')
axis equal; 
axis([0 4  0 1.7]);
set(gca,'xtick',0:1:4)
 xlabel('days in culture')
 ylabel('glucose(mg/ml)')
hold off

figure(3)

hold on
plot(time, dydt1(:,4),'-')
plot (time1 ,L_real,'v')
set(gca,'xtick',0:1:4)
%legend ('simulation', 'Data')
 axis equal; 
  axis([0 4  0 1.49]);
  xlabel('days in culture')
  ylabel('lactate(mg/ml)')
hold off










dydt1
end



 %  dydt (:,1)    ===>>>   p_0
 % P_0 = dydt(:,1)
% figure(1)
% hold on 
 %  plot ( time, p_0_real, '-o') % real date results
 %  plot (time , P_0, ' -*')     % simulation rsesults
%hold off 






%[p,fval,exitflag] = fminsearch(@leastcomp,p,[],td,hare,lynx);

function dydt = Mondoza(t,y,p) 


global tau_o tau_g p_star tau_og tau_go lambda L_star Kk_o K_L K_G Gamma n_star  m_star  G_star 


p_0 = y(1) ; %
p_g = y(2) ; %
G   = y(3) ; %
L   = y(4) ; %



Kar_L = 0.5 * (1+ tanh( Gamma *(y(4)-L_star))) ;

Kar_L_star = 1- Kar_L;

%%
%Kar G function 



G_min = 0.3;   % find it ! 
y(3);

 if (y(3) >G_min)
    Kar_G = 1;
 else
     Kar_G = 0 ; 

 end




%%

%5(a)
tmp1 =  (1/tau_o)*( 1- y(1) - y(2))* y(1)  + (1/tau_go).* Kar_L .* y(2)  - (1/tau_og).* Kar_L_star .*  Kar_G .* y(1);


tmp2 = (1/tau_g) .* (1- y(1) - y(2)) .* y(2)  - ( 1/tau_go).*Kar_L .* y(2) +(1/tau_og).* Kar_L_star .*  Kar_G .* y(1);


tmp3 =  - Kk_o  .*( y(3)./(y(3)+ lambda .* y(4) + n_star)).*y(1) -  K_G.*(y(3)./(y(3)+G_star)) .*y(2); 

%a= - Kk_o  .*( y(3)./(y(3)+ lambda .* y(4) + n_star)).*y(1)
%b = -  K_G.*(y(3)./(y(3)+G_star)) .*y(2)


tmp4 = - K_L  .*( y(4)./((y(3)./lambda)+  y(4) + m_star)).*y(1)     +   0.82*K_G.*(y(3)./(y(3)+G_star)) .*y(2); 
%c= - K_L  .*( y(4)./((y(3)./lambda)+  y(4) + m_star)).*y(1)




 dydt = [tmp1; tmp2; tmp3;  tmp4  ];

end
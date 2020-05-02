clearvars;close all;clc;
% Download the data from ref [1] and read them with the function getDataCOVID

time = datetime(2020,03,11):1:datetime(2020,05,1);
Confirmed = [178 240 301 379 466 568 716 867 1048 1263 ...
    1523 1822 2247 2637 3179 4048 4896 5662 6349 7388 8572 ...
    9965 11126 12565 13914 15193 16363 17625 19031 20169 21388 ...
    22416 23324 24372 25779 27098 28157 29425 30812 31530 ...
    33863 35842 37699 39561 41355 42626 43704 45200 46446 ...
    48829 50411 52238] ;
Recovered = [0 0 0 0 0 0 0 0 0 0 0 1 3 8 14 20 32 32 36 ...
    36 66 87 140 157 205 215 220 275 392 645 710 875 913 1051 ...
    1256 1259 1473 2560 2740 3053 3337 3345 3580 3602 3608 3608 ...
    3608 3636 3636 3636 3643 3643];
Deaths = [4 4 5 5 6 11 14 16 19 24 28 35 43 53 67 82 102 121 ...
    132 153 182 228 245 280 323 350 396 453 507 545 599 634 682 732 ...
    789 889 973 1057 1148 1180 1229 1322 1439 1533 1621 1698 1722 1788 ...
    1875 1956 2045 2136];

Npop= 39.56e6; % population


Recovered(1:2)=[];
Deaths(1:2)=[];
time(1:2)= [];
Confirmed(1:2)=[];

S0 = Npop-Confirmed(1);
E0 = 0;

alpha_hist = 1;

figure(1)
hold on

for WinID = 1:7
    
% time = datetime(2020,02,24):1:datetime(2020,03,17);
WinR=Recovered((WinID-1)*5+1:(WinID-1)*5+15);
WinD=Deaths((WinID-1)*5+1:(WinID-1)*5+15);
WinC=Confirmed((WinID-1)*5+1:(WinID-1)*5+15);
Wint=time((WinID-1)*5+1:(WinID-1)*5+15);



% Definition of the first estimates for the parameters
alpha_guess =0.07; % protection rate
beta_guess = 0.9; % Infection rate
gamma_guess = 1/2; 
lambda_guess = [0.02,0.05]; % recovery rate
kappa_guess = [0.02,0.05]; % death rate

guess = [alpha_guess,...
    beta_guess,...
    gamma_guess,...
    lambda_guess,...
    kappa_guess];


% Initial conditions
I0 = WinC(1)-WinR(1)-WinD(1); % Initial number of infectious cases. Unknown but unlikely to be zero.
R0 = WinR(1);
D0 = WinD(1);

% Parameter estimation with the lsqcurvefit function
[alpha1,beta1,gamma1,Lambda1,Kappa1] = ...
    fit_SEIRDP(WinC-WinR-WinD,WinR,WinD,Npop,S0,E0,time,Wint,guess);

dt = 0.1; % time step
time1 = datetime(Wint(1)):dt:datetime(Wint(1)+24);
N = numel(time1);
t = [round(datenum(Wint(1)-time(1))/dt):round(datenum(Wint(1)-time(1))/dt)+N-1].*dt;
[S,E,I,R,D,P,BRN] = SEIRDP(alpha1,beta1,gamma1,Lambda1,Kappa1,Npop,S0,E0,I0,R0,D0,t);

BRN0(WinID) = mean(BRN(1:140)*alpha_hist*(1-alpha1)^15);
alpha_hist = alpha_hist*(1-alpha1)^5;
TimeBRN(WinID) = datetime(Wint(1)+14);

Startpoint = round(datenum(6)/dt);
S0 = S(Startpoint);
E0 = E(Startpoint);

I_predict = I(round(datenum(15)/dt):round(datenum(19)/dt));
time_predict = datenum(Wint(1)+15-time(1)):dt:datenum(Wint(1)+19-time(1));
time_real = 0:numel(time)-1;
I_real = interp1(time_real,Confirmed - Recovered -Deaths,time_predict);
delta_I(WinID) = mean(I_predict-I_real)/ mean(I_real);
delta_time(WinID) = datetime(Wint(1)+19);

R_predict = R(round(datenum(15)/dt):round(datenum(19)/dt));
time_predict = datenum(Wint(1)+15-time(1)):dt:datenum(Wint(1)+19-time(1));
time_real = 0:numel(time)-1;
R_real = interp1(time_real,Recovered,time_predict);
delta_R(WinID) = mean(R_predict-R_real)/ mean(R_real);

D_predict = D(round(datenum(15)/dt):round(datenum(19)/dt));
time_predict = datenum(Wint(1)+15-time(1)):dt:datenum(Wint(1)+19-time(1));
time_real = 0:numel(time)-1;
D_real = interp1(time_real,Deaths,time_predict);
delta_D(WinID) = mean(D_predict-D_real)/ mean(D_real);


 semilogy(time1,I,'r--',time1,R,'b--',time1,D,'k--');
 
end

WinID = WinID +1; 
WinR=Recovered((WinID-1)*5+1:(WinID-1)*5+15);
WinD=Deaths((WinID-1)*5+1:(WinID-1)*5+15);
WinC=Confirmed((WinID-1)*5+1:(WinID-1)*5+15);
Wint=time((WinID-1)*5+1:(WinID-1)*5+15);

% Initial conditions
I0 = WinC(1)-WinR(1)-WinD(1); % Initial number of infectious cases. Unknown but unlikely to be zero.
R0 = WinR(1);
D0 = WinD(1);

% Parameter estimation with the lsqcurvefit function
[alpha1,beta1,gamma1,Lambda1,Kappa1] = ...
    fit_SEIRDP(WinC-WinR-WinD,WinR,WinD,Npop,S0,E0,time,Wint,guess);

dt = 0.1; % time step
time1 = datetime(Wint(1)):dt:datetime(Wint(1)+104);
N = numel(time1);
t = [round(datenum(Wint(1)-time(1))/dt):round(datenum(Wint(1)-time(1))/dt)+N-1].*dt;
[S,E,I,R,D,P,BRN] = SEIRDP(alpha1,beta1,gamma1,Lambda1,Kappa1,Npop,S0,E0,I0,R0,D0,t);

BRN0(WinID) = mean(BRN(140)*alpha_hist*(1-alpha1)^15);
alpha_hist = alpha_hist*(1-alpha1)^5;
TimeBRN(WinID) = datetime(Wint(1)+14);

Startpoint = round(datenum(6)/dt);
S0 = S(Startpoint);
E0 = E(Startpoint);

 Fig1=semilogy(time1,I,'r',time1,R,'b',time1,D,'k');

Fig2 = semilogy(time,Confirmed-Recovered-Deaths,'ro',time,Recovered,'bo',time,Deaths,'ko');



% ylim([0,1.1*Npop])
title('California,USA (5/1/2020)')
ylabel('Number of cases')
xlabel('time (days)')
leg = {'Infected (fitted)',...
        'Recovered (fitted)','Death (fitted)',...
        'Infected (reported)','Recovered (reported)','Death  (reported)'};

legend([Fig1(1),Fig1(2),Fig1(3),Fig2(1),Fig2(2),Fig2(3)],leg)
set(gcf,'color','w')
grid on
axis tight
% ylim([1,8e4])
 set(gca,'yscale','lin');
 

% Fig5=semilogy(datetime(time1(1))+t_seirdp,I_seirdp,'r--',datetime(time1(1))+t_seirdp,R_seirdp,'b--',datetime(time1(1))+t_seirdp,D_seirdp,'k--');

figure(2)
plot(TimeBRN,BRN0,'o-',TimeBRN,ones(size(TimeBRN,2)),'r')
title('California,USA (5/1/2020)')
ylabel('Effective Reproduction Number (R_t)')
xlabel('time (days)')



figure(3)
plot(delta_time,delta_I,'r',delta_time,delta_R,'b',delta_time,delta_D,'k');
title('California,USA (5/1/2020)')
ylabel('Relative errors of predictions')
xlabel('time (days)')
leg = {'Infected',...
        'Recovered','Death'};
legend(leg{:})
yticklabels({'-70%','-60%','-50%','-40%','-30%','-20%','-10%','0','10%','20%','30%'})

figure(4)
bar(time(2:(WinID-1)*5+15),Confirmed(2:(WinID-1)*5+15)-Confirmed(1:(WinID-1)*5+14));


time_predict = (datenum(Wint(1)+14-time(1)):1:datenum(Wint(1)+104-time(1)));
Confirmed_predict = interp1(t,I+R+D,time_predict);

hold on
time2 = datetime(time1(1)+15):1:datetime(time1(1)+104);
bar(time2,Confirmed_predict(2:91)-Confirmed_predict(1:90));



title('California,USA (5/1/2020)')
ylabel('Daily increase of the confirmed infected cases')
xlabel('time (days)')
leg = {'Real data',...
        'Projection by proposed method'};
legend(leg{:})

figure(5)

bar(time(2:(WinID-1)*5+15),Deaths(2:(WinID-1)*5+15)-Deaths(1:(WinID-1)*5+14));

time_predict = (datenum(Wint(1)+14-time(1)):1:datenum(Wint(1)+104-time(1)));
Deaths_predict = interp1(t,D,time_predict);

hold on
time2 = datetime(time1(1)+15):1:datetime(time1(1)+104);
bar(time2,Deaths_predict(2:91)-Deaths_predict(1:90));


title('California,USA (5/1/2020)')
ylabel('Daily increase of the death cases')
xlabel('time (days)')
leg = {'Real data',...
        'Projection by proposed method'};
legend(leg{:})

figure(6)
bar(time(1:(WinID-1)*5+15),Confirmed(1:(WinID-1)*5+15));

time_predict = (datenum(Wint(1)+15-time(1)):1:datenum(Wint(1)+104-time(1)));
Confirmed_predict = interp1(t,I+R+D,time_predict);

hold on
time2 = datetime(time1(1)+15):1:datetime(time1(1)+104);
bar(time2,Confirmed_predict);

title('California,USA (5/1/2020)')
ylabel('Number of the confirmed infected cases')
xlabel('time (days)')
leg = {'Real data',...
        'Projection by proposed method'};
legend(leg{:})

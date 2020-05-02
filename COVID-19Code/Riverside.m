clearvars;close all;clc;
% Download the data from ref [1] and read them with the function getDataCOVID

time = datetime(2020,03,19):1:datetime(2020,05,1);
Confirmed = [20 26 28 43 46 57 105 105 183 195 233 291 371 429 493 638 ...
    665 799 946 1016 1179 1280 1350 1431 1619 1751 1961 2105 2264 2457 ...
    2602 2638 2847 2960 3084 3218 3315 3409 3563 3643 3735 3942 4031 4079] ;
Recovered = [0 0 0 0 0 0 0 0 0 0 0 0 30 40 50 50 60 60 60 67 74 125 156 ...
    156 194 297 376 376 472 472 610 642 700 700 700 700 700 700 700 700 ...
    700 700 700 700];
Deaths = [3 4 5 6 6 6 8 8 8 8 8 9 13 13 14 15 18 19 25 28 32 33 39 ...
    41 41 50 50 54 59 69 74 75 85 93 99 100 112 117 118 141 141 143 ...
    149 156];

Npop= 2.471e6; % population


Recovered(1:4)=[];
Deaths(1:4)=[];
time(1:4)= [];
Confirmed(1:4)=[];

S0 = Npop-Confirmed(1);
E0 = 0;

alpha_hist = 1;

figure(1)
hold on

for WinID = 1:5
    
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
title('Riverside,CA,USA (5/1/2020)')
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
title('Riverside,CA,USA (5/1/2020)')
ylabel('Effective Reproduction Number (R_t)')
xlabel('time (days)')



figure(3)
plot(delta_time,delta_I,'r',delta_time,delta_R,'b',delta_time,delta_D,'k');
title('Riverside,CA,USA (5/1/2020)')
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



title('Riverside,CA,USA (5/1/2020)')
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


title('Riverside,CA,USA (5/1/2020)')
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

title('Riverside,CA,USA (5/1/2020)')
ylabel('Number of the confirmed infected cases')
xlabel('time (days)')
leg = {'Real data',...
        'Projection by proposed method'};
legend(leg{:})

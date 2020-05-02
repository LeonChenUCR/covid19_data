clearvars;close all;clc;
time = datetime(2020,03,11):1:datetime(2020,05,1);
Confirmed = [218 327 423 617 748 971 1409 2556 5699 8530 ...
    12323 16893 21727 26389 33038 39141 46124 53520 59454 67387 76051 ...
    84046 92743 103158 114173 123159 131829 140385 151078 161807 172972 ...
    181825 190295 196115 207168 218911 227625 235395 242580 248431 253325 ...
    258589 263744 271162 279314 287490 293696 298442 301489 306478 311379 ...
    315515] ;
Recovered = [14 14 14 14 14 14 14 137 215 317 467 607 757 1035 1485 2029 2689 3545 4185 ...
    5101 6276 7619 9185 10858 12790 14325 15615 17583 19950 22121 24180 25232 ...
    25798 26463 26917 27511 28089 28227 28884 29275 29593 29826 30018 31418 ...
    31962 32106 32344 32805 32904 33332 33678 34486];
Deaths = [0 0 0 2 6 10 16 17 34 56 75 148 184 269 366 461 596 834 ...
    966 1162 1398 1781 2326 2955 3549 4151 4758 5489 6268 7076 7845 8628 9385 ...
    10056 14612 15648 16130 17086 17627 18298 18822 19754 20255 20982 21430 ...
    21802 22365 22668 23134 23480 23796 24039];

Npop= 8.623e6; % population

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
title('New York,USA (5/1/2020)')
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
title('New York,USA (5/1/2020)')
ylabel('Effective Reproduction Number (R_t)')
xlabel('time (days)')



figure(3)
plot(delta_time,delta_I,'r',delta_time,delta_R,'b',delta_time,delta_D,'k');
title('New York,USA (5/1/2020)')
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



title('New York,USA (5/1/2020)')
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


title('New York,USA (5/1/2020)')
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

title('New York,USA (5/1/2020)')
ylabel('Number of the confirmed infected cases')
xlabel('time (days)')
leg = {'Real data',...
        'Projection by proposed method'};
legend(leg{:})

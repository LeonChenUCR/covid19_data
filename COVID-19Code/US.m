clearvars;close all;clc;
% Download the data from ref [1] and read them with the function getDataCOVID
time = datetime(2020,03,09):1:datetime(2020,05,1);
Confirmed = [583,959,1281,1663,2179,2727,3679,4644,6098,8758,14183,...
    19713,26878,35238,45069,55672,69658,86568,105361,125413,143442,...
    165224,190982,217540,247278,279869,313995,339785,370960,405673,...
    438916,473540,508455,538648,566186,593070,623517,654425,686222,...
    717835,746377,772215,800954,827250,856875,891270,928514,960132,...
    987215,1011692,1036930,1065771,1097242,1132292];
Recovered = [7,8,8,12,12,12,21,21,24,149,266,369,519,668,1099,1403,...
    2110,2810,3615,4688,5525,6693,8258,10043,12039,14480,16874,19180,...
    21490,24565,28197,32441,37426,40476,42485,44777,50087,52981,56386,...
    59941,66829,71003,73527,76048,93942,98133,101518,108394,109483,...
    114041,118410,126515,132733,141897];
Deaths = [22,28,36,40,47,54,73,100,123,154,210,277,350,513,646,834,...
    1167,1400,1768,2232,2533,3081,3742,4612,5587,6585,8113,9379,...
    10853,12886,14692,16599,18708,20599,22088,23639,26020,30922,...
    33308,36571,38098,40652,42237,45805,48032,50314,52334,54079,...
    55411,56724,59169,61605,63731,65597];

Npop= 3.272e8; % population

Recovered(1:4)=[];
Deaths(1:4)=[];
time(1:4)= [];
Confirmed(1:4)=[];


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
title('USA (5/1/2020)')
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
title('USA (5/1/2020)')
ylabel('Effective Reproduction Number (R_t)')
xlabel('time (days)')



figure(3)
plot(delta_time,delta_I,'r',delta_time,delta_R,'b',delta_time,delta_D,'k');
title('USA (5/1/2020)')
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



title('USA (5/1/2020)')
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


title('USA (5/1/2020)')
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

title('USA (5/1/2020)')
ylabel('Number of the confirmed infected cases')
xlabel('time (days)')
leg = {'Real data',...
        'Projection by proposed method'};
legend(leg{:})



function [S,E,I,R,D,P,BRN] = SEIRDP(alpha,beta,gamma,Lambda,Kappa,Npop,S0,E0,I0,R0,D0,t)
% [S,E,I,Q,R,D,P] = SEIQRDP(alpha,beta,gamma,delta,lambda,kappa,Npop,E0,I0,R0,D0,t)
% simulate the time-histories of an epidemic outbreak using a generalized
% SEIR model.
%
% Input
%
%   alpha: scalar [1x1]: fitted protection rate
%   beta: scalar [1x1]: fitted  infection rate
%   gamma: scalar [1x1]: fitted  Inverse of the average latent time
%   delta: scalar [1x1]: fitted  inverse of the average quarantine time
%   lambda: scalar [1x1]: fitted  cure rate
%   kappa: scalar [1x1]: fitted  mortality rate
%   Npop: scalar: Total population of the sample
%   E0: scalar [1x1]: Initial number of exposed cases
%   I0: scalar [1x1]: Initial number of infectious cases
%   Q0: scalar [1x1]: Initial number of quarantined cases
%   R0: scalar [1x1]: Initial number of recovered cases
%   D0: scalar [1x1]: Initial number of dead cases
%   t: vector [1xN] of time (double; it cannot be a datetime)
%
% Output
%   S: vector [1xN] of the target time-histories of the susceptible cases
%   E: vector [1xN] of the target time-histories of the exposed cases
%   I: vector [1xN] of the target time-histories of the infectious cases
%   Q: vector [1xN] of the target time-histories of the quarantinedcases
%   R: vector [1xN] of the target time-histories of the recovered cases
%   D: vector [1xN] of the target time-histories of the dead cases
%   P: vector [1xN] of the target time-histories of the insusceptible cases
%
% Author: E. Cheynet - UiB - last modified 16-03-2020
%
% see also fit_SEIQRDP.m

%% Initial conditions
N = numel(t);
Y = zeros(6,N);
Y(1,1) = S0;
Y(2,1) = E0;
Y(3,1) = I0;
Y(4,1) = R0;
Y(5,1) = D0;
Y(6,1) = Npop-S0-E0-R0-D0-I0;

if round(sum(Y(:,1))-Npop)~=0
    error('the sum must be zero because the total population (including the deads) is assumed constant');
end
%%
modelFun = @(Y,A,F) A*Y + F;
dt = median(diff(t));
% ODE resolution

lambda = Lambda(1)*(1-exp(-Lambda(2).*t)); % I use these functions for illustrative purpose only
kappa = Kappa(1)*exp(-Kappa(2).*t); % I use these functions for illustrative purpose only
BRN = beta./(lambda+kappa);

for ii=1:N-1
    A = getA(alpha,gamma,lambda(ii),kappa(ii));
    SI = Y(1,ii)*Y(3,ii);
    F = zeros(6,1);
    F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
    Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
end


S = Y(1,1:N);
E = Y(2,1:N);
I = Y(3,1:N);
R = Y(4,1:N);
D = Y(5,1:N);
P = Y(6,1:N);



    function [A] = getA(alpha,gamma,lambda,kappa)
        A = zeros(6);
        % S
        A(1,1) = -alpha;
        % E
        A(2,2) = -gamma;
        % I
        A(3,2:3) = [gamma,-kappa-lambda];
        % R
        A(4,3) = lambda;
        % D
        A(5,3) = kappa;
        % P
        A(6,1) = alpha;
    end
    function [Y] = RK4(Fun,Y,A,F,dt)
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*dt,A,F);
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end

end



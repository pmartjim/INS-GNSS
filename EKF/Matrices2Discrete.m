function [Phi, Gamma, Gammaw_W_Gammaw]= Matrices2Discrete( F, Gu_tilda, Gw, H, Qw, dT )
%MATRICES2DISCRETE Summary of this function goes here

% Detailed explanation goes here

sysc= ss(F, Gu_tilda, H, 0);
sysd= c2d(sysc, dT);
Phi= sysd.A;
Gamma= sysd.B;

% Methdo to obtain covariance matrix for dicrete system
C= [-F, Gw*Qw*Gw';
     zeros(15), F'];
 
EXP= expm(C*dT);
Gammaw_W_Gammaw= EXP(16:end,16:end)' * EXP(1:15,16:end);
% Gammaw_W_Gammaw= (Gw*dT) * (Qw/dT) * (Gw*dT)'; % simplified version


end


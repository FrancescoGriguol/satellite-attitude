%% Novembre 2023

% Gruppo 09: Marotti Dario, Majorana Andres Maria, Griguol Francesco, Scrascia Damiano


% Script di inizializzazione per la seconda esercitazione di CAS

clc
clear
close all


%% Prima esercitazione

% Dati numerici 

h=717;                                    % quota orbitale [km]
I_x=0.025;                                % momento di inerzia asse x [kg*m^2]
I_y=0.0255;                               % momento di inerzia asse y [kg*m^2]
I_z=0.007;                                % momento di inerzia asse z [kg*m^2]
I=[I_x; I_y; I_z];
mu_E=398600;                              % costante gravitazionale Terra [km^3/s^2]
R_E=6378.14;                              % raggio terrestre [km]
w_o=sqrt(mu_E/(h+R_E)^3);                 % velocità orbitale [rad/s]

% Dati coppie di disturbo

A_Fx=2e-8;                                % ampiezza oscillazione coppia in direzione x [N/m]
A_Fy=1e-8;                                % ampiezza oscillazione coppia in direzione y [N/m]
A_Fz=4e-8;                                % ampiezza oscillazione coppia in direzione z [N/m]

f_Fx=5e-2;                                % pulsazione coppia in direzione x [rad/s]
f_Fy=3e-2;                                % pulsazione coppia in direzione y [rad/s]
f_Fz=2e-2;                                % pulsazione coppia in direzione z [rad/s]

ph_Fx=pi/2;                               % fase coppia in direzione x [rad]
ph_Fy=pi/2;                               % fase coppia in direzione y [rad]
ph_Fz=pi/2;                               % fase coppia in direzione z [rad]


w_ic=[0; -w_o; 0];
qi=[1; 0; 0; 0];

%% Seconda Esercitazione

% Specifiche di sistema

e_rp = 1e-2;                              % [rad]
tr = 60;                                  % [s]
m_phase = 60;                             % [deg]
w_A = log(10)/tr;                         % [s^-1]


alpha_phi = 1/(w_o^2*(I_z-I_y));        
alpha_psi = 1/(w_o^2*(I_y-I_x));

% Funzioni di trasferimento della dinamica

G_phi = tf(-alpha_phi,[-alpha_phi*I_x 0 1]);
G_theta = tf(1,[I_y 0 0]);
G_psi = tf(alpha_psi,[alpha_psi*I_z 0 1]);

% Guadagni del controllore

K_phi = 1/abs(alpha_phi)*(1/e_rp -1);
K_theta = 6e-6;
K_psi = 1/abs(alpha_psi)*(1/e_rp -1);

% Verifica della stabilità

if alpha_phi > 0 && K_phi > 1/abs(alpha_phi)
    disp('Stabilità verificata per rollio');
elseif alpha_phi < 0 && K_phi > 0
    disp('Stabilità verificata per rollio')
end

if K_theta > 0
    disp("Stabilità verificata per beccheggio")
else
    disp("Stabilità non verificata per beccheggio")
end

if alpha_psi > 0 && K_psi > 1/abs(alpha_psi)
    disp('Stabilità verificata per imbardata');
elseif alpha_psi < 0 && K_psi > 0
    disp('Stabilità verificata per imbardata')
end

% Funzioni di trasferimento dei controllori ottenute dal ControlSystemDesigner

load("C_phi.mat")
load("C_theta.mat")
load("C_psi.mat")

% Definizione del controllore con Control System Designer
% Funzioni di trafserimento della catena aperta

CG_phi = C_phi*G_phi;
CG_theta = C_theta*G_theta;
CG_psi = C_psi*G_psi;

%Tracciamento dei diagrammi di Bode

figure(1)
hold on
grid on
bode(CG_phi,'r')
title('TF catena aperta roll')

figure(2)
hold on
grid on 
bode(CG_theta,'b')
title('TF catena aperta pitch')

figure(3)
hold on
grid on 
bode(CG_psi,'g')
title('TF catena aperta yaw')

% Funzioni di trasferimento del feedback

W_phi = feedback(CG_phi,1);
W_theta = feedback(CG_theta,1);
W_psi = feedback(CG_psi,1);

% Tracciamento della risposta al gradino
step_value = 20;                           % [deg]
opt = stepDataOptions;
opt.StepAmplitude = step_value;

[phi,tphi] = step(W_phi,opt);
[theta,ttheta] = step(W_theta,opt);
[psi,tpsi] = step(W_psi,opt);

trphi = interp1(phi,tphi,0.9*step_value,'pchip','extrap');
trtheta = interp1(theta,ttheta,0.9*step_value,'pchip','extrap');
trpsi = interp1(psi,tpsi,0.9*step_value,'pchip','extrap');

t = 0:0.1:60 ;
figure(4)
subplot(3,1,1);
hold on; box on; grid on;
step(20*W_phi,t,'r')
yline(step_value,'k--');
yline(0.9*step_value,'k--');
xline(trphi,'k--',['t_r = ',num2str(trphi),' s']);
title('\phi(t)');
xlabel('[s]'); ylabel('[deg]'); ylim([0 25]);

subplot(3,1,2);
hold on; box on; grid on;
step(20*W_theta,t,'-b')
yline(step_value,'k--');
yline(0.9*step_value,'k--');
xline(trtheta,'k--',['t_r = ',num2str(trtheta),' s']);
title('\theta(t)');
xlabel('[s]'); ylabel('[deg]');ylim([0 25]);

subplot(3,1,3);
hold on; box on; grid on;
step(20*W_psi,t,'-g')
yline(step_value,'k--');
yline(0.9*step_value,'k--');
xline(trpsi,'k--',['t_r = ',num2str(trpsi),' s']);
title('\psi(t)');
xlabel('[s]'); ylabel('[deg]');ylim([0 25]);






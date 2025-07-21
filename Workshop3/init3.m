% ESERCITAZIONE 3
% -------------------------------------------------------------------------
%
% Gruppo 9:
% Marotti Dario-Majorana Andres Maria-Griguol Francesco-Scrascia Damiano
% 
% File di inizializzazione del modello Simulink

clc
clear
close all

% CARATTERISTICHE DELL'ORBITA E CONDIZIONI INIZIALI

mu = 398600;                   % [km^3/s^2] costante gravitazionale Terra
h = 717;                       % [km]       altitudine orbita  
R_E = 6378.14;                 % [km]       raggio equatoriale
r_o = h+R_E;                   % [km]       raggio orbitale
w_o = sqrt(mu/(h+R_E)^3);      % [rad/s]    velocità orbitale 
w_ic=[0; -w_o; 0];             % [rad/s]    velocità orbitale 
i = 0;                         % [°]        inclinazione orbitale

% Momento di dipolo residuo del satellite
d_sp = [0.001; 0.003; 0.005];

% Condizioni iniziali dinamica di Eulero
y_0= [0; -w_o; 0; 0; 0];

% Quaternione iniziale
q_0= [1; 0; 0; 0];

% CARATTERISTICHE DEI PANNELLI SOLARI - ENDUROSAT 3U SOLAR PANEL
% (I pannelli sono modellati come una massa singola concentrata 
%  all'estremità di un'asta dotata di rigidezza e smorzamento)

mp = 0.146;                    % [kg] per pannello solare
Lp = 0.2;                      % [m]  lunghezza caratteristica
kt = 0.002;                    % []   tra 0.002 e 0.02
kp = kt/(Lp^2);                % []   coeff. di rigidezza
eps = 0.005;                   % []   smorzamento critico
cp = 2*eps*sqrt(kp*mp);        % []   coeff. di smorzamento

% MOMENTI DI INERZIA 

% Momenti di inerzia del satellite Cubesat 3U
I_x = 0.025;                     % [kg*m^2] momento di inerzia asse x
I_y = 0.0255;                    % [kg*m^2] momento di inerzia asse y 
I_z = 0.007;                     % [kg*m^2] momento di inerzia asse z 

% Momenti di inerzia totali comprendenti i pannelli solari (nominalmente
% nel piano x-y del sistema body)

I_y = I_y + 2*mp*Lp^2;         % [kg*m^2] momento di inerzia asse y
I_z = I_z + 2*mp*Lp^2;         % [kg*m^2] momento di inerzia asse z 

I = [I_x; I_y; I_z];           % vettore dei momenti di inerzia



% RUOTE DI REAZIONE 
% Modello ruota di reazione: CubeWheel S azienda Cubespace
% Modello motore DC brushless: Faulhaber, series 0620...B

Hm = 1.7e-3;                   % [Nms]    momento angolare massimo
wmax = 8000*pi/30;             % [rad/s]  velocità massima
Tmax = 0.28e-3;                % [Nm]     coppia massima
Iw = Hm/wmax;                  % [Kg*m^2] momento d'inerzia della ruota
Km= 1.09e-3;                   % [Nm/A]   costante di coppia
Kv = Km;                       % [Nm/V]   coef di tensione del motore
Rm = 8.8;                      % [Ohm]    resistenza elettrica del motore
B = 1.02e-9/60;                % [Nm/s]   coefficiente di attrito viscoso

% Guadagno di retroazione sulla corrente K >> Kv*Km/Iw
% (si moltiplica per un fattore 20)
Kprp = Km*Kv/Iw*20;           

% NUOVE FUNZIONI DI TRASFERIMENTO DEL SISTEMA

% Parametri delle funzioni di trasferimento del sistema
alpha_phi = 1/(4*w_o^2*(I_z-I_y));
alpha_theta = 1/(3*w_o^2*(I_x-I_z));
alpha_psi = 1/(w_o^2*(I_y-I_x));

% Funzioni di trasferimento del sistema modificato per l'effetto di
% gradiente gravitazionale

G_phi = tf(-alpha_phi,[-alpha_phi*I_x 0 1]);
G_theta = tf(alpha_theta,[alpha_theta*I_y 0 1]);
G_psi = tf(alpha_psi,[alpha_psi*I_z 0 1]);

% Condizioni di stabilità sui guadagni del controllore
e_rp = 1e-2;         % [rad]
tr = 60;             % [s]
m_phase = 60;        % [deg]
w_A = log(10)/tr;    % [1/s]

K_phi = 1/abs(alpha_phi)*(1/e_rp -1);
K_theta = 1/abs(alpha_theta)*(1/e_rp -1);
K_psi = 1/abs(alpha_psi)*(1/e_rp -1);

if alpha_phi > 0 && K_phi > 1/abs(alpha_phi)
    disp('Stabilità verificata per rollio');
elseif alpha_phi < 0 && K_phi > 0
    disp('Stabilità verificata per rollio')
end

if alpha_theta > 0 && K_theta > 1/abs(alpha_theta)
    disp("Stabilità verificata per beccheggio")
else
    disp("Stabilità non verificata per beccheggio")
end

if alpha_psi > 0 && K_psi > 1/abs(alpha_psi)
    disp('Stabilità verificata per imbardata');
elseif alpha_psi < 0 && K_psi > 0
    disp('Stabilità verificata per imbardata')
end

% Funzioni di trasferimento del controllore

load("C_phi.mat")
load("C_theta.mat")
load("C_psi.mat")

% % REIEZIONE DEI DISTRUBI A BASSA FREQUENZA
% 
% d0 = 5e-8;               % [Nm]    ampiezza del disturbo
% wtilde = 1.05*w_o;       % [rad/s] frequenza del disturbo
% ed = e_rp;               % errore max dovuto all'ingresso di disturbo
% 





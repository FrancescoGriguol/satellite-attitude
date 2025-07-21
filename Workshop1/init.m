%Script di inizializzazione per la prima esercitazione di CAS
clc
clear
close all


%%
%Dati numerici

h=717;                                       %quota orbitale [km]
I_x=0.028;                                 %momento di inerzia asse x [kg*m^2]
I_y=0.024;                                 %momento di inerzia asse y [kg*m^2]
I_z=0.007;                                 %momento di inerzia asse z [kg*m^2]
I=[I_x; I_y; I_z];
mu_E=398600;                          %costante gravitazionale Terra [km^3/s^2]
R_E=6378.14;                           % raggio terrestre [km]
w_o=sqrt(mu_E/(h+R_E));       %velocità orbitale [km/s]

%Dati coppie di disturbo

A_Fx=2e-8;                               %ampiezza oscillazione coppia in direzione x [N/m]
A_Fy=1e-8;                               %ampiezza oscillazione coppia in direzione y [N/m]
A_Fz=4e-8;                               %ampiezza oscillazione coppia in direzione z [N/m]

f_Fx=5e-3;                                %pulsazione coppia in direzione x [rad/s]
f_Fy=9e-3;                                %pulsazione coppia in direzione y [rad/s]
f_Fz=2e-2;                                %pulsazione coppia in direzione z [rad/s]

ph_Fx=pi/2;                              % fase coppia in direzione x [rad]
ph_Fy=pi/2;                              % fase coppia in direzione x [rad]
ph_Fz=pi/2;                              % fase coppia in direzione x [rad]


w_ic=[0; -w_o; 0];
qi=[1; 0; 0; 0];

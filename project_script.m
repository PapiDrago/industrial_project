close all;
clc;
clear;
% Fatto fino al punto1
Ld = 0.05;     % [H]
Lq = 0.011;    % [H]
Rs = 0.62;     % [Ohm]
we = 785;      % [rad/s]
vdc = 650;     % [V]

% === MATRICE A ===
A = [-Rs/Ld,     we*Lq/Ld;
     -we*Ld/Lq,  -Rs/Lq];
olc_eigenvalues = eig(A);

% === MATRICE B ===
B = [1/Ld,     0;
     0,        1/Lq];

% === USCITA E MATRICE D ===
C = eye(2);
D = zeros(2);

x0_1 = 0; %-15
x0_2 = 0;%+15
x0 = [x0_1;x0_2];


u1_max = vdc/sqrt(3);
u2_max = vdc/sqrt(3);
u1_min = -vdc/sqrt(3);
u2_min = -vdc/sqrt(3);
u_sat_max = [u1_max; u2_max];
u_sat_min = [u1_min; u2_min];

x1_max = 15;
x2_max = 15;
x1_min = -15;
x2_min = -15;
x_sat_max = [x1_max; x2_max];
x_sat_min = [x1_min; x2_min];

% === SISTEMA CONTINUO ===
sys_c = ss(A, B, C, D);

% Reachability control
Co = ctrb(A, B);
rank_C = rank(Co);

% Observability control
Ob = obsv(A, eye(2));
rank_O = rank(Ob);

% === DISCRETIZZAZIONE ===
Ts = 1/8000;  % Periodo di campionamento [s]
sys_d = c2d(sys_c, Ts);
[Ad, Bd, Cd, Dd] = ssdata(sys_d);
old_eigenvalues = eig(Ad);
% Reachability control
Cod = ctrb(Ad, Bd);
rank_Cd = rank(Cod);

% Observability control
Obd = obsv(Ad, eye(2));
rank_Od = rank(Obd);

xref = [6.9; 8.9];
n = length(xref);
nu = size(B, 2);

%% === LQR CONTINUO ===
Q = 1*eye(2);  % Penalizza errori sugli stati
R = 1*eye(2);  % Penalizza l'uso eccessivo del controllo

u_bar_c = -inv(B)*A*xref;
[Kc,Sc,Pc] = lqr(A, B, Q, R);
eig_lqr = eig(A-B*Kc);
%% === LQR DISCRETO ===
Qd = 1000*eye(2);
Rd = 1*eye(2); 
%Qd2 = 1000*eye(2);
%Rd2= 1*eye(2); 
N=5;
N2=5;
N3=100;

u_bar_d = Bd\ (eye(size(Ad))-Ad)*xref;

[Kd,Sd,Pd] = dlqr(Ad, Bd, Qd, Rd);
%[Kd2,Sd2,Pd2] = dlqr(Ad, Bd, Qd2, Rd);
eig_lqr = eig(Ad-Bd*Kd);
%open('regulator_LQ_project')
%sim('regulator_LQ_project')
%keyboard

S0=zeros(2);

%% Point 5

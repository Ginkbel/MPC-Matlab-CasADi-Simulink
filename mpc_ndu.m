%% NDU Seminar
% Jurica Vuckovic; JMBAG: 0035219927

clear all
clc
close all


A = [-0.03125 0.625 0.0025;-1 0 0.08;0.16 -3.2 -0.0128];
B = [0;0;0.1];
C = [1 0 0];
D = 0;

[num,den] = ss2tf(A,B,C,D);

sys = tf(num,den);
pole(sys);

t_val = 0:1:249;

sys_ss = ss(A,B,C,D);
y1 = step(sys_ss, t_val);

ts = timeseries(1*ones(250,1));

out = sim('ndu');

figure(1)
step(t_val,sys_ss)
title('Odziv modela prostora stanja i Simulink modela sustava')
xlabel('t [s]'); ylabel('\omega [rad/s]') 
hold on
grid on
plot(out.w2{1}.Values, 'r')
legend('Prostor stanja','Simulink model')
hold off

% Ekvivalentan zapis; ali drugaciji su brojevi, odzivi su isti
[a,b,c,d] = linmod('ndu');

sys_simulink = ss(a,b,c,d);
y2 = step(sys_simulink, t_val);

tocnost = sum(y1 - y2) % Provjera koliko su blizu odzivi funkcija  (5.1792e-13 rezultat -- jednaki odzivi)

% out = sim('ndu');

figure(2)
plot(ts, '--r','LineWidth',1.5)
title('Odziv Simulink sustava na jediničnu pobudu')
xlabel('t [s]'); ylabel('\omega [rad/s]') 
hold on
grid on
plot(out.w2{1}.Values, 'b','LineWidth',1.5)
legend('m_1_,_r_e_f','\omega_2')
hold off

% Diskretizacija 

Td = 0.5; % vrijeme uzorkovanja
[Ad,Bd,Cd,Dd] = c2dm(A,B,C,D,Td,'zoh');

sys_ssd = c2d(sys_ss,Td);

figure(3)
step(sys_ssd,10, 'b')
hold on
grid on
step(sys_ss,10, 'r')
legend('Diskretizirani model prostora stanja', 'Kontinuirani model prostora stanja')

% MPC - formulacija

addpath('D:\Faks 3. god\R2021a\bin\casadi-windows-matlabR2016a-v3.5.5') % change path to your directory containing Casadi
import casadi.*

Td = 0.5;
N = 150; % predikcijski horizont

% Ogranicenja
m1_max = 200; m1_min = -m1_max; % maksimalni i minimalni ulazni moment
alfa_max = 0.2; alfa_min = - alfa_max; % maksimalni i minimalni kut torzije (x2 varijabla stanja)

% Definicija varijabli stanja modela
omega2 = SX.sym('w2'); del_alfa = SX.sym('del_alfa'); omega1 = SX.sym('w1');
states = [omega2;del_alfa;omega1]; n_states = length(states);

m1 = SX.sym('m1');
control = m1; n_controls = length(control);

jdzba = Ad*states+Bd*control; % vrijednost varijabli stanja u koraku (k+1)
f = Function('f',{states,control},{jdzba}); % linearna funkcija f_k+1(x_k,u_k)

U = SX.sym('U',n_controls, N); % upravljacke varijable kroz cijeli vremenski horizont
P = SX.sym('P',n_states + 1);
% parametri; pocetna stanja varijabli stanja sustava i referentna brzina
% vrtnje tereta (x_1 varijabla stanja)

X = SX.sym('X',n_states,(N+1));
% matrica koja sadrži stanja sustava preko cijelog horizonta N (dim. 3xN)

X(:,1) = P(1:3); % početno stanje sustava
for k = 1:N
    st = X(:,k);  con = U(:,k);
    f_value = f(st,con);
    X(:,k+1) = f_value;
end

ff=Function('ff',{U,P},{X}); % za spremanje optimalnih predvidenih vrijednosti

obj = 0; % fja cilja (objective function)
g = [];  % vektor ogranicenja vrijednosti

q = 100; % tezinski faktor razlike brzina (u sumi)
r = 0.1; % tezinski faktor ulazne velicine
qn = 1; % tezinski faktor terminal cost-a


% Funkcija cilja - suma
for k = 1:N
    st = X(1,k);  con = U(:,k);
    obj = obj + q *((st - P(4))^2) + r*(con^2); % Suma
end

% Dodaje se jos clan terminal cost na funkciju cilja
obj = obj + qn*((X(1,N+1)-P(4))^2);

% Koristiti ce se za ogranicenje |del_alfa| <= 0.2rad
for k = 1:N+1   
    g = [g ; X(2,k)];   % del_alfa ogranicenje
end

% definiranje optimizacijskih varijabli
OPT_variables = reshape(U,N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


args = struct;

% Ogranicenje |del_alfa| <= 0.2rad 
args.ubg = alfa_max; args.lbg = alfa_min;

% Ogranicenje |m1| <= 200Nm
args.ubx = m1_max; args.lbx = m1_min;

% Glavna petlja simulacije MPC problema
%--------------------------------------

t0 = 0;
x0 = [0 ; 0 ; 0.0];    % Pocetni uvjeti stanja sustava [omega2;del_alfa;omega1]
xs = 5; % Referentna brzina vrtnje vretena (omega2) [rad/s]

xx(:,1) = x0; % xx matrica sadrzi vrijednosti stanja sustava za svaki korak optimizacije
t(1) = t0; % vektor vremena

u0 = zeros(N,1);  % vektor svih optimizacijskih varijabli; dim. Nx1

% sim_tim = 20; % Maksimalno vrijeme simulacije

% Pocetak MPC problema
iter = 0;
xx1 = [];
u_cl= [];

main_loop = tic; % Pocetak brojanja vremena

while (norm((x0(1)-xs),2) > 1e-2 && iter < 2500)
    args.p   = [x0;xs]; % pocetne vrijednosti stanja sustava i referentna brzina vrtnje
    args.x0 = u0; % pocetne vrijednosti optimizacijskih varijabli

    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    
    u = full(sol.x);
    ff_value = ff(u',args.p); % racunanje predikcije stanja na temelju optimalnih koraka upravljackog vektora u
    
    xx1(:,1:3,iter+1)= full(ff_value)'; % spremanje prva 3 predvidena stanja u 3D matricu
    
    u_cl = [u_cl ; u(1,:)]; % Spremanje svih upravljackih signala u vektor (samo prva vrijednost optimalnog rjesenja upravljackog svakog koraka)
    
    t(iter+1) = t0; % Spremanje u vektor vremena
    
    [t0, x0, u0] = shift(Td, t0, x0, u,f); % inicijalizacija sljedeceg koraka
    
    xx(:,iter+2) = x0; % Spremanje stanja sustava u matricu stanja  
    
    iter = iter + 1;
end


main_loop_time = toc(main_loop) % zaustavlja brojanje vremena i ispisuje koliko je ukupno trajala petlja

ss_error = norm((x0(1)-xs),2)
average_computation_mpc_step = main_loop_time/iter
fprintf('Broj iteracija MPC petlje: %d \n' ,iter)

figure(4)

plot(t,xs*(ones(length(t),1)),'--b','LineWidth',1.5)
ylim([0 1.05*xs])
hold on
stairs(t,xx(1,2:size(xx,2))','k','LineWidth',1.75)
grid on
tit = title('Odziv \omega_2'); tit.FontSize = 16;
leg = legend('\omega_2_R','\omega_2');set(leg, 'FontSize',13);set(leg, 'Location','southeast')
x_lab = xlabel('t [s]'); y_lab = ylabel('\omega [rad/s]');
set(x_lab, 'FontSize',15), set(y_lab, 'FontSize',15)

figure(5)

subplot(2,1,1)
stairs(t,u_cl','r','LineWidth',1.5)
tit = title('Upravljačka varijabla m_1 i \Delta\alpha'); tit.FontSize = 16;
grid on
leg = legend('m_1 [Nm]'); set(leg, 'FontSize',11);
x_lab = xlabel('t [s]'); y_lab = ylabel('m_1 [Nm]');
set(x_lab, 'FontSize',13), set(y_lab, 'FontSize',13)

subplot(2,1,2)
stairs(t,xx(2,2:size(xx,2))', 'g', 'LineWidth',1.75)
ylim([-inf 0.22])
grid on
leg = legend('\Delta\alpha [rad]'); set(leg, 'FontSize',11);
x_lab = xlabel('t [s]'); y_lab = ylabel('\Delta\alpha [rad]');
set(x_lab, 'FontSize',13); set(y_lab, 'FontSize',13)

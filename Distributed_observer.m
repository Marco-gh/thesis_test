clc;clear;close all;rng('default')
%% DEF SIS
% Matrice instabile
% A = [1 2 -2 -2;1.5 -0.5 1.5 0;1.5 0.5 -2.5 -2;2 2 -2 -3];
% Moti oscillatori e stabili
% A = [-0.1   -1.5    0      0;
%       1.5   -0.1    0      0;
%       0      0     -0.05  -0.8;
%       0      0      0.8   -0.05];
% Moti oscillatori e instabili
A = [ 0.10  -1.50   0      0;
    1.50   0.10   0      0;
    0      0      0.05  -0.80;
    0      0      0.80   0.05];

C{1} = [1 0 0 1];
C{2} = [0 1 1 0];
C{3} = [0 0 1 0];
C{4} = [0 0 0 1];

[V,D,W] = eig(A);
n = size(A,1);
%% Test osservabilità globale
H = vertcat(C{:});
assert(rank(obsv(A,H))==n)
%% Def rete
n_agents = 4;
Lgraph = [2 -1 0 -1 
     0 1 -1 0
     -1 -1 2 0
     -1 0 0 1];
Adj = diag(diag(Lgraph))-Lgraph;
%% Matrici di osservabilità
rankO = cell(1,n_agents);
O = cell(1,n_agents);
for i=1:n_agents
    O{i} = obsv(A,C{i});
    rankO{i} = rank(O{i});
end
%% DECOMPOSIZIONE DI OSSERVABILITA'
Abar = cell(1,n_agents);
Bbar = cell(1,n_agents);
Cbar = cell(1,n_agents);
T = cell(1,n_agents);
for i=1:n_agents
    [Abar{i},Bbar{i},Cbar{i},T{i}]=detectabilty_decomposition(O{i},A,zeros(n,1),C{i});
end
%% Calcolo guadagni da "A Simple Approach to Distributed Observer Design for Linear Systems"
% Punti 1-2 pg. 5
r = null(Lgraph');
r = n_agents*r/sum(r);
R = diag(r);

% Grafo fortemente connesso
assert(rank(Lgraph) == n_agents-1)
assert(all(r>0))

Lhat = R*Lgraph+Lgraph'*R;

g = cell(1,n_agents);
for i = 1:n_agents
    rankO_i = rankO{i};

    g{i} = [eye(rankO_i), zeros(rankO_i,n-rankO_i);
            zeros(n-rankO_i,rankO_i), zeros(n-rankO_i)];
end
G = blkdiag(g{:});
T_d = blkdiag(T{:});

% Punto 3
tol = 1e-10;
step = 1e-5;
eps = 0;
Acheck = T_d' * kron(Lhat,eye(n)) * T_d + G;
while min(eig(Acheck - (eps + step)*eye(n*n_agents))) > tol
    eps = eps + step;
end

% Punto 4
alpha = 1e-4; % Tasso di decadimento
gamma = cell(1,n_agents); % Coefficienti di accoppiamento globale
for i=1:n_agents
    gamma{i}=0;
    v = rankO{i};
    nu = n - v;
    
    % Per nodo completamente osservabile
    if nu == 0
        gamma{i}=0;
        continue
    end

    A_iu = Abar{i}(v+1:end,v+1:end);
    A_ir = Abar{i}(v+1:end,1:v);
    
    gamma{i} = step;
    
    while true
        d = gamma{i}*eps/r(i) - 2*alpha;

        % Controllo denominatore positivo
        if d <= 0
            gamma{i} = gamma{i} + step;
            continue
        end
        
        Acheck = A_iu+A_iu'-d*eye(nu)+1/d*(A_ir*A_ir');
        if max(eig(Acheck))<0
            break
        end
        gamma{i} = gamma{i} + step;
    end
end
% q positivo con un margine del 10 percernto
gamma=max(max(cell2mat(gamma)),1.1*max(2*alpha*r/eps));

% Punto 5
L_io = cell(1,n_agents);
A_io = cell(1,n_agents);
C_io = cell(1,n_agents);
for i=1:n_agents
    v = rankO{i};
    if v>0
        A_io{i} = Abar{i}(1:v,1:v);
        C_io{i} = Cbar{i}(:,1:v);

        poles_desired = -0.5*(1:v);

        L_io{i} = place(A_io{i}', C_io{i}', poles_desired)';
        % Per verifica
        eig_check = eig(A_io{i} - L_io{i}*C_io{i});
    else
        L_io{i} = zeros(0,size(Cbar{i},1));
    end
end

% Punto 6-7 - Individuare Q di Lyapunov e guadagno locale
P_io = cell(1,n_agents);
L = cell(1,n_agents);
M = cell(1,n_agents);
for i=1:n_agents
    v = rankO{i};
    nu=n-v;

    F = A_io{i} - L_io{i}*C_io{i} + alpha*eye(v);
    q = gamma*eps/r(i) - 2*alpha;
    P_io{i} = lyap(F', q*eye(v));
    assert(min(eig(P_io{i}))>0)

    p_i = size(C_io{i},1); % Numero di uscite locali

    L{i} = T{i}*[L_io{i};zeros(nu,p_i)];
    M{i} = T{i}*[inv(P_io{i}) zeros(v,nu);zeros(nu,v) eye(nu)]*T{i}';
end
%% Def stato esteso, parametri simulazione comuni
xd=@(t,x) A*x;
x0 = randn(n,1);
xhat0 = cell(1,n_agents);
% Vettore di stato con stato originale + tutte le stime (n+n*n_agents):
z0 = x0;
for i = 1:n_agents
    xhat0{i} = randn(n,1)*2.5;
    % xhat0{i} = zeros(n,1);

    z0 = [z0;xhat0{i}];
end

%% Parametri simulazione
tmax = 150;
tau = 0.5;
Phi = expm(A*tau); % Usato in più punti
%% Sim sensa ritardi
% [t,z]=ode15s(@(t,z) distributedObserver_ODE(z,A,C,L,M,gamma,Adj,n_agents,n),[0 tmax],z0);
%% Sim con delay senza predizione
in_history = @(t)z0;
sol = dde23(@(t,z,Z)distributedObserver_DDE(z,Z,A,C,L,M,gamma,Adj,n_agents,n),tau,in_history,[0 tmax]);
t_np_obs = sol.x';
z_np_obs = sol.y';
%% Sim con delay e predizione
in_history = @(t)z0;
sol = dde23(@(t,z,Z)distributedObserver_DDE_predictor(z,Z,A,C,Phi,L,M,gamma,Adj,n_agents,n),tau,in_history,[0 tmax]);
t_pred_obs = sol.x';
z_pred_obs = sol.y';
%% Sim dinamiche globali errore con delay
% Senza compensazione
Lambda_cell = cell(1,n_agents);
for i = 1:n_agents
    Lambda_cell{i} = A - L{i}*C{i};
end
Lambda = blkdiag(Lambda_cell{:});
Mglob  = blkdiag(M{:});
A0_nc = Lambda;
A1_nc = -gamma*Mglob*kron(Lgraph,eye(n));
e0 = zeros(n*n_agents,1);
for i = 1:n_agents
    idx = (i-1)*n + (1:n);
    e0(idx) = xhat0{i} - x0;
end
sol_nc_err = dde23(@(t,e,E) A0_nc*e + A1_nc*E(:,1),tau,e0,[0 tmax]);
t_nc_err = sol_nc_err.x';
e_np_err = sol_nc_err.y';
% Con compensazione
Dgraph = diag(sum(Adj,2));
A0_pred = Lambda-gamma*Mglob*kron(Dgraph,eye(n));
A1_pred = gamma*Mglob*kron(Adj,Phi);
sol_err = dde23(@(t,e,E) A0_pred*e + A1_pred*E(:,1),tau,e0,[0 tmax]);
t_pred_err = sol_err.x';
e_pred_err = sol_err.y';
%% Plot dinamiche agente
% agent_to_show = 3;
% idx = n + (agent_to_show-1)*n + (1:n);
% % Senza compensazione
% figure;subplot(2,1,1);hold on;grid on
% plot(t_np_obs,z_np_obs(:,idx))
% title('Without compensation')
% hold off
% % Con compensazione
% subplot(2,1,2);hold on;grid on
% plot(t_pred_obs,z_pred_obs(:,idx))
% title('With predictor')
% hold off
%% Plot dinamiche errori
figure;hold on;grid on;
% Senza compensazione
subplot(2,1,1);hold on;grid on
for i = 1:n_agents
    idx = (i-1)*n + (1:n);
    plot(t_nc_err,vecnorm(e_np_err(:,idx),2,2));
end
xlabel('t')
ylabel('||e_i||_2')
title('Without compensation')
legend('e_1','e_2','e_3','e_4')
hold off
% Con compensazione
subplot(2,1,2);hold on;grid on
for i = 1:n_agents
    idx = (i-1)*n + (1:n);
    plot(t_pred_err,vecnorm(e_pred_err(:,idx),2,2));
end
xlabel('t')
ylabel('||e_i||_2')
title('With predictor')
legend('e_1','e_2','e_3','e_4')
hold off
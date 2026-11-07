clc;clear;close all;rng('default')
%% Parametri simulazione
tmax = 500;
tau = 5;
alpha = 1e-3; % Tasso di decadimento
%% DEF SIS
% 3 modi marginalmente stabili
A = [0 1 0 0 0 0 
     -1 0 0 0 0 0
     0 0 0 2 0 0
     0 0 -2 0 0 0
     0 0 0 0 0 3
     0 0 0 0 -3 0];
C{1}=[1 0 0 0 0 0]; % Vede modi 1,2
C{2}=[0 0 1 0 0 0]; % Vede modo 2
C{3}=[0 0 0 0 0 1]; % Vede modo 3

[~,eigA,~] = eig(A);
n = size(A,1);
%% Stato aumentato [x; u_obs]
A_aug = [A eye(n);zeros(n) zeros(n)];
exp_aug = expm(A_aug*tau);
Phi = exp_aug(1:n,1:n);
Gamma = exp_aug(1:n,n+1:end);
%% Def rete
n_agents = length(C);
% Grafo fortemente connesso diretto
% Lgraph=[3 -1 -1 -1
%        -1 1 0 0
%        -1 0 1 0
%        -1 0 0 1];
% Grafo NON diretto, L simmetrica
Lgraph = [2 -1 -1
          -1 1 0
          -1 0 1];
Adj = diag(diag(Lgraph))-Lgraph;
%% Test osservabilità globale
H = vertcat(C{:});
rank(H)
assert(rank(obsv(A,H))==n)
%% Matrici di osservabilità
rankO = cell(1,n_agents);
O = cell(1,n_agents);
for i=1:n_agents
    O{i} = obsv(A,C{i});
    rankO{i} = rank(O{i});
end
%% DECOMPOSIZIONE DI OSSERVABILITA
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

% Controllo per grafo fortemente connesso
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
Acheck = T_d'*kron(Lhat,eye(n))*T_d+G;
while min(eig(Acheck-(eps+step)*eye(n*n_agents))) > tol
    eps = eps + step;
end

% Punto 4
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
        d = gamma{i}*eps/r(i)-2*alpha;

        % Controllo denominatore positivo
        if d <= 0
            gamma{i} = gamma{i}+step;
            continue
        end
        
        Acheck = A_iu+A_iu'-d*eye(nu)+1/d*(A_ir*A_ir');
        if max(eig(Acheck))<0
            break
        end
        gamma{i} = gamma{i}+step;
    end
end
% q positivo con un margine del 10 percernto
gamma=max(max(cell2mat(gamma)),1.0*max(2*alpha*r/eps));

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

        L_io{i} = place(A_io{i}',C_io{i}',poles_desired)';
        % Per verifica
        eig_check = eig(A_io{i}-L_io{i}*C_io{i});
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

    F = A_io{i}-L_io{i}*C_io{i}+alpha*eye(v);
    q = gamma/r(i) - 2*alpha;
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
%% Sim sensa ritardi
% [t,z]=ode15s(@(t,z) distributedObserver_ODE(z,A,C,L,M,gamma,Adj,n_agents,n),[0 tmax],z0);
%% Sim con delay senza predizione
in_history = @(t)z0;
sol = dde23(@(t,z,Z)distributedObserver_DDE(z,Z,A,C,L,M,gamma,Adj,n_agents,n),tau,in_history,[0 tmax]);
t_np_obs = sol.x';
z_np_obs = sol.y';
%% Sim con delay e predizione 1 (solo propagazione libera)
in_history = @(t)z0;
sol = dde23(@(t,z,Z)distributedObserver_DDE_predictor_1(z,Z,A,C,Phi,L,M,gamma,Adj,n_agents,n),tau,in_history,[0 tmax]);
t_pred_obs_1 = sol.x';
z_pred_obs_1 = sol.y';
%% Sim con delay e predizione 2 (con correzione costante)
sol = dde23(@(t,z,Z)distributedObserver_DDE_predictor_2(z,Z,A,C,Phi,Gamma,L,M,gamma,Adj,n_agents,n),tau,in_history,[0 tmax]);
t_pred_obs_2 = sol.x';
z_pred_obs_2 = sol.y';
%% Sim dinamiche globali errore con delay
Lambda_cell = cell(1,n_agents);
for i = 1:n_agents
    Lambda_cell{i} = A-L{i}*C{i};
end
Lambda = blkdiag(Lambda_cell{:});
Mglob  = blkdiag(M{:});
e0 = zeros(n*n_agents,1);
for i = 1:n_agents
    idx = (i-1)*n+(1:n);
    e0(idx) = xhat0{i}-x0;
end
% Senza predizione
A0_np = Lambda;
A1_np = -gamma*Mglob*kron(Lgraph,eye(n));
sol_nc_err = dde23(@(t,e,E) A0_np*e+A1_np*E(:,1),tau,e0,[0 tmax]);
t_nc_err = sol_nc_err.x';
e_np_err = sol_nc_err.y';
% Con predizione 1
Dgraph = diag(sum(Adj,2));
A0_pred = Lambda-gamma*Mglob*kron(Dgraph,eye(n));
A1_pred = gamma*Mglob*kron(Adj,Phi);
sol_err = dde23(@(t,e,E) A0_pred*e+A1_pred*E(:,1),tau,e0,[0 tmax]);
t_pred_err = sol_err.x';
e_pred_err_1 = sol_err.y';
% Con predizione 2
Dgraph = diag(sum(Adj,2));
Bpred_cell = cell(1,n_agents);
for j = 1:n_agents
    Bpred_cell{j} = Phi - Gamma*L{j}*C{j};
end
Bpred = blkdiag(Bpred_cell{:});
% Dinamica globale errore
A0_pred2 = Lambda-gamma*Mglob*kron(Dgraph,eye(n));
A1_pred2 = gamma*Mglob*kron(Adj,eye(n))*Bpred;
% Simulazione DDE
sol_err2 = dde23(@(t,e,E) A0_pred2*e + A1_pred2*E(:,1),tau,e0,[0 tmax]);
t_pred_err_2 = sol_err2.x';
e_pred_err_2 = sol_err2.y';
%% Plot dinamiche agente
agent_to_show = 1;
idx = n+(agent_to_show-1)*n+(1:n);
figure;hold on;
% Senza predizione
subplot(3,1,1);grid on
plot(t_np_obs,z_np_obs(:,[1 2 3 4]))
plot(t_np_obs,z_np_obs(:,idx),'--')
title('Without prediction')
legend('$x_1$','$x_2$','$x_3','$x_4','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$','Interpreter','latex')
hold off
% Con predizione 1
subplot(3,1,2);hold on;grid on
plot(t_pred_obs_1,z_pred_obs_1(:,[1 2 3 4]))
plot(t_pred_obs_1,z_pred_obs_1(:,idx),'--')
title('With predictor')
legend('$x_1$','$x_2$','$x_3','$x_4','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$','Interpreter','latex')
hold off
% Con predizione 2
subplot(3,1,3);hold on;grid on
plot(t_pred_obs_2,z_pred_obs_2(:,[1 2 3 4]))
plot(t_pred_obs_2,z_pred_obs_2(:,idx),'--')
title('With predictor')
legend('$x_1$','$x_2$','$x_3','$x_4','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$','Interpreter','latex')
hold off
%% Plot dinamiche errori
figure;hold on;grid on;
% Senza predizione
subplot(3,1,1);hold on;grid on
for i = 1:n_agents
    idx = (i-1)*n+(1:n);
    plot(t_nc_err,vecnorm(e_np_err(:,idx),2,2));
end
xlabel('t')
ylabel('||e_i||_2')
title('Without prediction')
legend('$e_1$','$e_2$','$e_3$','$e_4$','Interpreter','latex')
hold off
% Con predizione 1
subplot(3,1,2);hold on;grid on
for i = 1:n_agents
    idx = (i-1)*n+(1:n);
    plot(t_pred_err,vecnorm(e_pred_err_1(:,idx),2,2));
end
xlabel('t')
ylabel('||e_i||_2')
title('With prediction')
legend('$e_1$','$e_2$','$e_3$','$e_4$','Interpreter','latex')
hold off
% Con predizione 3
subplot(3,1,3);hold on;grid on
for i = 1:n_agents
    idx = (i-1)*n+(1:n);
    plot(t_pred_err_2,vecnorm(e_pred_err_2(:,idx),2,2));
end
xlabel('t')
ylabel('||e_i||_2')
title('With prediction')
legend('$e_1$','$e_2$','$e_3$','$e_4$','Interpreter','latex')
hold off
%% Comparazione errori finali
for i = 1:n_agents
    idx = (i-1)*n + (1:n);
    fprintf('Pred1 e%d final = %.4e\n',i,norm(e_pred_err_1(end,idx)));
    fprintf('Pred2 e%d final = %.4e\n',i,norm(e_pred_err_2(end,idx)));
end
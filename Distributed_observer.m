clc;clear;close all
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
% Serve a garantire q positivo nel calcolo sotto con un margine del 10
% percernto
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
    M{i} = T{i}*[inv(P_io{i}) zeros(v,nu);zeros(nu,v) eye(nu)]*T{i}'
end

%% SIM PLANT E NODI
tmax = 150;

xd=@(t,x) A*x;
x0 = randn(n,1);
xhat0 = cell(1,n);
% Vettore di stato con stato originale + tutte le stime (n+n*n_agents):
z0 = x0;
for i = 1:n_agents
    xhat0{i} = randn(n,1)*2.5;
    % xhat0{i} = zeros(n,1);

    z0 = [z0;xhat0{i}];
end

% Con le tolleranze di ode45 esplodeva
[t,z]=ode15s(@(t,z) distributedObserver_ODE(z,A,C,L,M,gamma,Adj,n_agents,n),[0 tmax],z0);
%% RES
x=z(:,1:n);
xhat=cell(1,n_agents);
for i=1:n_agents
    idx=n+(i-1)*n+(1:n);
    xhat{i}=z(:,idx);
end
%% PLOT STATO E USCITA
figure;hold on;grid on
for i=1:n_agents
    e=xhat{i}-x;
    % plot(t,e)
    plot(t,vecnorm(e,2,2)) % Norma 2 fatte per righe (una per x_i)
end
legend('e1','e2','e3','e4')
xlabel('t')
ylabel('||xhat_i-x||'); hold off
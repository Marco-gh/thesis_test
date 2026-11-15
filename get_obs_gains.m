function [L,M,gamma] = get_obs_gains(Lgraph,n,n_agents,O,T,Abar,Cbar,rankO,alpha)
%   Calcolo guadagni da "A Simple Approach to Distributed Observer Design for Linear Systems"
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
end
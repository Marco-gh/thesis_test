function dz = distributedObserver_DDE_predictor_2(z,Z,A,C,Phi,Gamma,L,M,gamma,Adj,n_agents,n)
    z_tau = Z(:,1);
    x = z(1:n);
    x_tau = z_tau(1:n);
    dz = zeros(size(z));
    dz(1:n) = A*x;

    %% Predizione per ogni nodo j
    % Quella affetta da ritardo e che verrà poi effettivamente passata, per
    % semplicità tutto calcolato qua dentro
    xhat_pred = cell(1,n_agents);
    for j = 1:n_agents
        idx_j = n+(j-1)*n+(1:n);
        xhat_j_tau = z_tau(idx_j);
        y_j_tau = C{j}*x_tau;
        innovation_j_tau = L{j}*(y_j_tau-C{j}*xhat_j_tau);
        xhat_pred{j} = Phi*xhat_j_tau+ Gamma*innovation_j_tau; % Approssimazione che verrà poi condivisa sulla rete
    end
    %% Dinamica di ogni observer i
    for i = 1:n_agents
        idx_i = n+(i-1)*n+(1:n);
        xhat_i = z(idx_i);
        y_i = C{i}*x;
        coupling = zeros(n,1);
        for j = 1:n_agents
            if Adj(i,j) ~= 0
                coupling = coupling+Adj(i,j)*(xhat_pred{j}-xhat_i);
            end
        end
        dz(idx_i) = A*xhat_i+L{i}*(y_i-C{i}*xhat_i)+gamma*M{i}*coupling;
    end
end
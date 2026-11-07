function dz = distributedObserver_DDE_predictor_1(z,Z,A,C,Phi,L,M,gamma,Adj,n_agents,n)
    z_tau = Z(:,1);
    x = z(1:n);

    dx = A*x;
    dz = zeros(size(z));
    dz(1:n) = dx;

    for i = 1:n_agents
        idx_i = n + (i-1)*n+(1:n);
        xhat_i = z(idx_i);
        y_i = C{i}*x;
        coupling = zeros(n,1);
        for j = 1:n_agents
            if Adj(i,j) ~= 0
                idx_j = n+(j-1)*n+(1:n);
                xhat_j_pred = Phi*z_tau(idx_j); % Predizione al tempo corrente
                coupling = coupling+Adj(i,j)*(xhat_j_pred-xhat_i);
            end
        end
        dxhat_i = A*xhat_i+L{i}*(y_i-C{i}*xhat_i)+gamma*M{i}*coupling;
        dz(idx_i) = dxhat_i;
    end
end
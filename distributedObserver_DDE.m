function dz = distributedObserver_DDE(z,Z,A,C,L,M,gamma,Adj,n_agents,n)
    z_tau = Z(:,1);
    x = z(1:n); % Stato reale corrente

    dx = A*x;
    dz = zeros(size(z));
    dz(1:n) = dx;

    for i = 1:n_agents
        % xhat_i(t)
        idx_i = n + (i-1)*n+(1:n);
        xhat_i = z(idx_i);
        xhat_i_tau = z_tau(idx_i);
        y_i = C{i}*x;
        coupling = zeros(n,1);
        for j = 1:n_agents
            if Adj(i,j) ~= 0
                idx_j = n+(j-1)*n+(1:n);
                xhat_j_tau = z_tau(idx_j);
                coupling = coupling+Adj(i,j)*(xhat_j_tau-xhat_i_tau);
            end
        end
        dxhat_i = A*xhat_i+L{i}*(y_i-C{i}*xhat_i)+gamma*M{i}*coupling;
        dz(idx_i) = dxhat_i;
    end
end
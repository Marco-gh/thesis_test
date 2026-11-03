function dz = distributedObserver_ODE(z,A,C,Lobs,M,gamma,Adj,n_agents,n)
    x = z(1:n);
    xhat = cell(1,n_agents);

    for i = 1:n_agents
        idx = n + (i-1)*n + (1:n);
        xhat{i} = z(idx);
    end

    % SISTEMA
    dx = A*x;

    % OBS
    dxhat = cell(1,n_agents);

    for i = 1:n_agents
        yi = C{i}*x;
        innovation = yi-C{i}*xhat{i};
        coupling = zeros(n,1);

        for j = 1:n_agents
            coupling = coupling+Adj(i,j)*(xhat{j}-xhat{i});
        end

        dxhat{i} = A*xhat{i}+Lobs{i}*innovation+gamma*M{i}*coupling;
    end

    %% Ricostruzione vettore derivata
    dz = dx;

    for i = 1:n_agents
        dz = [dz; dxhat{i}];
    end
end
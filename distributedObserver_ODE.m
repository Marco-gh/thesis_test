function dz = distributedObserver_ODE(z,A,C,L,M,gamma,Adj,N,n)
    x = z(1:n);
    xhat = cell(1,N);

    for i = 1:N
        idx = n + (i-1)*n + (1:n);
        xhat{i} = z(idx);
    end

    % SISTEMA
    dx = A*x;

    % OBS
    dxhat = cell(1,N);

    for i = 1:N
        yi = C{i}*x;
        innovation = yi-C{i}*xhat{i};
        coupling = zeros(n,1);

        for j = 1:N
            coupling = coupling+Adj(i,j)*(xhat{j}-xhat{i});
        end

        dxhat{i} = A*xhat{i}+L{i}*innovation+gamma*M{i}*coupling;
    end

    %% Ricostruzione vettore derivata
    dz = dx;

    for i = 1:N
        dz = [dz; dxhat{i}];
    end
end
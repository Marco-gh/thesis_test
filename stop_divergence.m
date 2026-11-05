function [value,isterminal,direction]=stop_divergence(t,z,Z,N,type,n,emax)
assert(type==1 || type==2)

errmax=0;

if type==1
    x=z(1:n);
    for i=1:N
        idx=n+(i-1)*n+(1:n);
        errmax=max(errmax,norm(z(idx)-x));
    end
else
    for i=1:N
        idx=(i-1)*n+(1:n);
        errmax=max(errmax,norm(z(idx)));
    end
end

value=emax-errmax;
isterminal=1;
direction=-1;
end
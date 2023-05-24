function nbV = newtonBasisAtX(x,centers,type,par,scale)
% Determines the Newton Basis at a given x for an RKHS
% The RKHS is defined by
% centers = a set of centers
% par = the RBF parameters and 
% scale = scaling factor 
N = length(centers);
V = zeros(N,N);
knlV = real(kermat(centers,x,type,par,scale));

Z = real(kermat(centers,centers,type,par,scale));
z = diag(Z);

v_x = zeros(N,1);
P = sqrt(real(kermat(centers(1),centers(1),type,par,scale)));

v_x(1) = knlV(1)/P;


if N > 1
    V(:,1) = real(kermat(centers(1),centers,type,par,scale))/P;
    w = V(:,1).^2;
    NN = V(2,:)*v_x;
    P = sqrt(z(2)-w(2));
    
    for jj = 1:N-1
        e = zeros(N,1);
        e(jj+1) = 1;
        Vm = V(:,1:jj);
        u = real(kermat(centers,centers(jj+1),type,par,scale))-Vm*Vm'*e;
        u_x = knlV(jj+1) - NN;
        v_x(jj+1) = u_x/P; % Newton Basis evaluated at x
        V(:,jj+1) = u/P; % Newton Basis evaluated at the centers
        w = V(:,jj+1).^2 + w; % Updating sum of Newton basis squared at the centers
        if jj ~=N-1
            NN = V(jj+2,:)*v_x;
            P = sqrt(z(jj+2)-w(jj+2));
        end
    end
end

nbV = v_x;

end
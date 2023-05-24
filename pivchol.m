function [R,P] = pivchol(A)
n = size(A,1); 
m = size(A,2);
if n ~= m
    error('Matrix is not symmetric')
end

R = zeros(n,m);

P = diag(ones(n,1));

for kk = 1:n
    B = A(kk:n,kk:n);
    maximum = max(diag(B));
    if kk ~= n
        [~,ll]=find(diag(A)==maximum);
        if ll > kk
            ll = ll(1);
            A(:,[ll kk]) = A(:,[kk ll]);
            R(:,[ll kk]) = R(:,[kk ll]);
            A([ll kk],:) = A([kk ll],:);
            P(:,[ll kk]) = P(:,[kk ll]);
        end
    end
    R(kk,kk) = sqrt(A(kk,kk));
    R(kk,kk+1:n) = 1/R(kk,kk)*A(kk,kk+1:n);
    A(kk+1:n,kk+1:n) = A(kk+1:n,kk+1:n)- R(kk,kk+1:n)'*R(kk,kk+1:n);
    
end
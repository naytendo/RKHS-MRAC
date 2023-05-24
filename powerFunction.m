function pow = powerFunction(type,x,kv,K,scale,par)

% Normalizaed
% pow = sqrt(kermat(x,x,type,par,scale) - kv'*(K\kv))/sqrt(kermat(x,x,type,par,scale));

% Not normalized
pow = sqrt(kermat(x,x,type,par,scale) - kv'*(K\kv));

%% Probably Outdated
% pow = sqrt(kernel(type,x,x,scale,par) - kv'*(K\kv));
end

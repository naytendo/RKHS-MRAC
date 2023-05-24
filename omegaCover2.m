function omega = omegaCover2(Ri,Ro,h)
[Z1,Z2]=meshgrid(-Ro-h-eps:h:Ro,-Ro-h-eps:h:Ro);
Z=[Z1(:),Z2(:)];
normZ = sqrt(sum(Z.^2,2));
Zp = Z(normZ < Ro,:);
normZp = sqrt(sum(Zp.^2,2));
omega = Zp(normZp > Ri,:);



end
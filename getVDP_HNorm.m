function [hnorm_fs,condK,Ns] = getVDP_HNorm(type,RBFpar,scale,mu,Ri,Ro)
h = [1;0.5;0.1;0.05];
Ns = zeros(length(h),1);
hnorm_fs = zeros(length(h),1);
condK = zeros(length(h),1);
for hh = 1:length(h)
    centers = omegaCover2(Ri,Ro,h(hh));
    Kc = real(kermat(centers,centers,type,RBFpar,scale));
    yEval = (1+mu)*centers(:,2)-mu*centers(:,1).^2.*centers(:,2);   
    theta =Kc\yEval;
    Ns(hh) = length(centers);
    %     norm_fs(index) = sqrt(coefs(1:Ns(nn))'*ur(1:Ns(nn),indexes)*coefs(1:Ns(nn)));
    %    
    hnorm_fs(hh) = sqrt(theta'*Kc*theta);
    condK(hh) = cond(Kc);
end



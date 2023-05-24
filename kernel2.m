function kern = kernel2(x,y,type,par,scale)
% Calculate value of kernel given the type of kernel k(.,.),two elements
% and parameters.
% Input:
%       type    -   string, type of kernel 
%       x,y     -   column vectors
%       beta    -   column vectors, parameters of kernel (indications to be added)
% Output:
%       KK      -   k(x,y)
r = norm(x-y);
s = r^2/2/scale^2;
% Examine the inputs (to be added)
switch type
    case 'g'
        Dxy = r^2;
        kern = 1/((sqrt(2*pi)*scale)^length(x))*exp(-1/2/scale^2 * Dxy);
    case 'wendland10'
        if 1-r/scale < 0
            kern = 0;
        else
            kern = (1-r/scale);
        end
     case 'w1'
        if 1-r/scale < 0
            kern = 0;
        else
            kern = (1-r/scale)^3*(3*r/scale+1);
        end
     case 'wendland12'
        if 1-r/scale < 0
            kern = 0;
        else
            kern = (1-r/scale)^5*(8*(r/scale)^2 + 5*r/(scale) + 1);
        end
    case 'w3'
        if 1-r/scale < 0
            kern = 0;
        else
            kern = (1-r/scale)^4*(4*r/scale+1);
        end
     case 'wendland32'
        if 1-r/scale < 0
            kern = 0;
        else
            kern = (1-r/scale)^6*(35*(r/scale)^2 + 18*r/scale +3);
        end
     case 'wendland33'
        if 1-r/scale < 0
            kern = 0;
        else
            kern = (1-r/scale)^8*(32*r^3/(scale^3) + 25*r^2/(scale^2) + 8*r/scale + 1);
        end
    case 'imq'
        kern = 1/(5^2+r^2)^scale;
    case 'matern32'
        kern = (1 + sqrt(3)*r/scale)*exp(-sqrt(3)*r/scale);
    case 'matern52'
        kern = (1 + sqrt(5)*r/scale + 5/3*r^2/(scale^2))*exp(-sqrt(5)*r/scale);
    case 'matern72'
        kern = (1+ sqrt(7)*r + 7/5*r^2/(scale^2)+7/15*sqrt(7)*r^3/(scale^3))*exp(-sqrt(7)*r/scale);
    case 'ms'
        nu = par;
        kern = 2^(1-nu)/gamma(nu)*(sqrt(2*nu)*r/scale+eps)^nu*besselk(nu,sqrt(2*nu)*r/scale+eps);
%         if r > 1e-17
%             sigma = 1;
%             kern = sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*r/hyper)^nu*besselk(nu,sqrt(2*nu)*r/hyper);
%         else
%             kern = 1;
%         end

    otherwise 
   
        error('Kernel type is not supported!');
end
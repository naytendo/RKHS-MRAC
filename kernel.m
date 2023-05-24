function kern = kernel(type,x,y,hyper,nu)
% Calculate value of kernel given the type of kernel k(.,.),two elements
% and parameters.
% Input:
%       type    -   string, type of kernel 
%       x,y     -   column vectors
%       beta    -   column vectors, parameters of kernel (indications to be added)
% Output:
%       KK      -   k(x,y)
r = norm(x-y)+eps;
% Examine the inputs (to be added)
switch type
    case {'g','exp'}
        Dxy = r^2;
        kern = 1/((sqrt(2*pi)*hyper)^length(x))*exp(-1/2/hyper^2 * Dxy);
    case 'wendland10'
        if 1-r/hyper < 0
            kern = 0;
        else
            kern = (1-r/hyper);
        end
     case {'wendland11','w1'}
        if 1-r/hyper < 0
            kern = 0;
        else
            kern = (1-r/hyper)^3*(3*r/hyper+1);
        end
     case 'wendland12'
        if 1-r/hyper < 0
            kern = 0;
        else
            kern = (1-r/hyper)^5*(8*(r/hyper)^2 + 5*r/(hyper) + 1);
        end
    case {'wendland31','w3'}
        if 1-r/hyper < 0
            kern = 0;
        else
            kern = (1-r/hyper)^4*(4*r/hyper+1);
        end
     case 'wendland32'
        if 1-r/hyper < 0
            kern = 0;
        else
            kern = (1-r/hyper)^6*(35*(r/hyper)^2 + 18*r/hyper +3);
        end
     case 'wendland33'
        if 1-r/hyper < 0
            kern = 0;
        else
            kern = (1-r/hyper)^8*(32*r^3/(hyper^3) + 25*r^2/(hyper^2) + 8*r/hyper + 1);
        end
    case 'imq'
        kern = 1/(5^2+r^2)^hyper;
    case 'matern32'
        kern = (1 + sqrt(3)*r/hyper)*exp(-sqrt(3)*r/hyper);
    case 'matern52'
        kern = (1 + sqrt(5)*r/hyper + 5/3*r^2/(hyper^2))*exp(-sqrt(5)*r/hyper);
    case 'matern72'
        kern = (1+ sqrt(7)*r + 7/5*r^2/(hyper^2)+7/15*sqrt(7)*r^3/(hyper^3))*exp(-sqrt(7)*r/hyper);
    case {'matern','ms'}
        sigma = 1;
        kern = sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*r/hyper)^nu*besselk(nu,sqrt(2*nu)*r/hyper);
%         if r > 1e-17
%             sigma = 1;
%             kern = sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*r/hyper)^nu*besselk(nu,sqrt(2*nu)*r/hyper);
%         else
%             kern = 1;
%         end

    otherwise 
   
        error('Kernel type is not supported!');
end
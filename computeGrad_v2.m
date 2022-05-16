% x: variables
% x = [x1;y1;x2;y2;...;xi;yi;...;xn;yn];
% and it can be simplified as x = [z1;z2;...zi;...;zn], where zi = [xi;yi]
% p: surface sample points
% pn: normal vector of point p
% R: radius of station
function [g1,g2,g3,g4] = computeGrad_v2(x,p,pn,R,k,t)
   n =length(x);
   m = length(p);
   m_ = length(pn);
   if(n<4)
       error('size of x = %d: x must greater or equal to 4 (at least 2 staitons).',n);
   end
   
   if(mod(n,2))
       error('size of x = %d: x must be a even number.',n);
   end

   if(m~=m_)
       error('size of p = %d, size of pn = %d : size must match.',m,m_);
   end

   if(mod(m,2))
       error('size of m = %d: m must be a even number.',n);
   end
   
   z = @(x,idx) [x((idx-1)*2+1);x((idx-1)*2+2)];
   % gradient of objective function
   %g1 = 2*(diag([ones(1,2),2*ones(1,n-4),ones(1,2)])+diag(-ones(1,n-2),2)+diag(-ones(1,n-2),-2))*x; 
   g1 = 2*(2*eye(n)+circshift(-eye(n),-1)+circshift(-eye(n),1))*x;
   
   % gradient for covering constraints
   g2 = zeros(n,m/2);
   const1 = R^2;
   for j = 1:(m/2)
       p_j = z(p,j);
       %pn_j = z(pn,j);
       for i = 1:(n/2)
           x_i = z(x,i);
           u = x_i-p_j;
           %h = 1/(1+exp(-k*u'*pn_j/norm(u)+const1-u'*u));
           %g2((i-1)*2+1:(i-1)*2+2,j) = -h*(1-h)*k*(pn_j/norm(u)-u'*pn_j*u/norm(u)^3-2*u);
           h = 1/(1+exp(-k*(const1-u'*u)));
           g2((i-1)*2+1:(i-1)*2+2,j) = h*(1-h)*k*(2*u);
            %g2((i-1)*2+1:(i-1)*2+2,j) = (exp(-k*(const1 - u'*u))+1)*2*u;
       end
   end

   % gradient for radius constraints
   g3 = zeros(n,n/2);
   for i = 1:((n/2)-1)
       index = (i-1)*2;
       x1 = z(x,i);
       x2 = z(x,i+1);
       g3(index+1:index+2,i) = 2*x1-2*x2;
   end
   x1 = z(x,n/2);
   x2 = z(x,1);
   g3(end-1:end,end) = 2*x1-2*x2;


   % gradient for boundary constraints
   b = 10; % it is 5, added 1 offset
   g4 = zeros(n,n/2);
   for i=1:2:(n/2)
       x_i = x(i);
       y_i = x(i+1);
    if(abs(x_i)>abs(y_i))
        g4(i,i) = -sign(x_i)/(b-abs(x_i));
        %g4(i,i) = -sign(x_i);
    elseif(abs(x_i)<abs(y_i))
        g4(i+1,i) = -sign(y_i)/(b-abs(y_i));
        %g4(i+1,i) = -sign(y_i);
    else
        g4(i:i+1,i) = -[sign(x_i)/(5-abs(x_i));sign(y_i)/(5-abs(y_i))]/2;
        %g4(i:i+1,i) = -[sign(x_i);sign(y_i)]/2;
    end
   end



end
    
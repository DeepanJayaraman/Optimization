function [xopt,fopt,niter,gnorm,dx] = Conjugate_Gradient_Branin_hoo_function(varargin)

% f - objective function
% g - gradient
% al - alpha function

if nargin==0
    % define starting point
    x0 = [2 1]';
elseif nargin==1
    % if a single input argument is provided, it is a user-defined starting
    % point.
    x0 = varargin{1};
else
    error('Incorrect number of input arguments.')
end

tol = 1e-6;% termination tolerance

maxiter = 10000;% maximum number of allowed iterations

dxmin = 1e-6;% minimum allowed perturbation

gnorm = inf; x = x0; niter = 0; dx = inf;% initialize gradient norm, optimization vector, iteration counter, perturbation

% f = @(x1,x2) (x2-(5.1/(4*pi^2))*x1^2+5*x1/pi-6)^2+10*(1-1/(8*pi))*cos(x1)+10;% define the objective function:
f =@(x1,x2) 12.096*x1^2+21.504*x2^2-1.7321*x1-x2;
f2 = @(x) f(x(1),x(2));% redefine objective function syntax for use with optimization:


figure(1); clf; ezcontour(f,[-10 10 -10 10 ]); axis equal; hold on% plot objective function contours for visualization:



d= - double(G(x(1),x(2)));

% Conjugate gradient algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = double(G(x(1),x(2)));
    %     d= - grad(x);
    gnorm = norm(g);
    
    % step size calculation - GOLDEN SECTION METHOD
    
    a=0;                            % start of interval
    b=2;                            % end of interval
    epsilon=0.000001;               % accuracy value
    iter= 50;                       % maximum number of iterations
    tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
    k=0;                            % number of iterations
    
    
    alpha1=a+(1-tau)*(b-a);             % computing alpha values
    alpha2=a+tau*(b-a);
    
    f_alpha1=F_alpha(alpha1,x(1),x(2),d(1),d(2));         % computing values in alpha points
    f_alpha2=F_alpha(alpha2,x(1),x(2),d(1),d(2));
    
    while ((abs(b-a)>epsilon) && (k<iter))
        k=k+1;
        if(f_alpha1<f_alpha2)
            b=alpha2;
            alpha2=alpha1;
            alpha1=a+(1-tau)*(b-a);
            
    f_alpha1=F_alpha(alpha1,x(1),x(2),d(1),d(2));         % computing values in alpha points
    f_alpha2=F_alpha(alpha2,x(1),x(2),d(1),d(2));
        else
            a=alpha1;
            alpha1=alpha2;
            alpha2=a+tau*(b-a);
            
    f_alpha1=F_alpha(alpha1,x(1),x(2),d(1),d(2));         % computing values in alpha points
    f_alpha2=F_alpha(alpha2,x(1),x(2),d(1),d(2));
        end
        
        k=k+1;
    end
    % chooses minimum alpha
    if(f_alpha1<f_alpha2)
        alpha = alpha1;
    else
        alpha = alpha2;
    end
    
    % take step:
    xnew = x + alpha*d;
    
    gnew = G(xnew(1),xnew(2));
    gnorm_new = norm(gnew);
    
    d_new = -gnew + (gnorm_new/gnorm)^2*d;
    
    d= d_new;
    
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    
    % plot current point
    plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    title('Conjugate gradient method')
    xlabel('\bf\it{x_1}')
    ylabel('\bf\it{x_2}')
    refresh
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew
    
    
    
    
end
xopt = x
fopt = f2(xopt)
niter = niter - 1
hold on;
plot(xopt(1),xopt(2),'*r')
end
% define the gradient of the objective
% function g = grad(x)
% g = [(2 * (x(2)^1 - (5.1 / (4*pi^2))*(x(1))^2 + (5/pi)*(x(1))^1 - 6))*(-2*(5.1 / (4*pi^2))*(x(1))+(5/pi)) - 10*(1-(1 / (8*pi)))*sin(x(1)^1);
%     (2 * (x(2)^1 - (5.1 / (4*pi^2))*(x(1))^2 + (5/pi)*(x(1))^1 - 6))];
% end
% function al = alphaFunc(alpha,x,d)
% %((x(2)-alpha*g(2))-(5.1/(4*pi^2))*(x(1)-alpha*g(1))^2+5*(x(1)-alpha*g(1))/pi-6)^2+10*(1-1/(8*pi))*cos((x(1)-alpha*g(1)))+10
% al = ((x(2)+alpha*d(2))-(5.1/(4*pi^2))*(x(1)+alpha*d(1))^2+5*(x(1)+alpha*d(1))/pi-6)^2+10*(1-1/(8*pi))*cos((x(1)+alpha*d(1)))+10;
% end

% function f = F(x1,x2)
% f = 12.096*x1^2+21.504*x2^2-1.7321*x1-x2;
% end
function g = G(x1,x2)
g = [ (3024*x1)/125 - 17321/10000;(5376*x2)/125 - 1];
end
function f_al = F_alpha(alpha, x1, x2, d1, d2) 
f_al =(1512*(x1 + alpha*d1)^2)/125 - x2 - (17321*alpha*d1)/10000 - alpha*d2 - (17321*x1)/10000 + (2688*(x2 + alpha*d2)^2)/125;
end
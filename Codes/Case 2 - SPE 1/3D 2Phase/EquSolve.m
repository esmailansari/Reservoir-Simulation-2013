function [x,Landa_W,Landa_O] = EquSolve(f,x0,TolX,MaxIter,varargin)
%% Defining the initial condition 
h = 1e-4; TolFun = 10^-6; EPS = 1e-6;
fx = feval(f,x0,varargin{:});
Nf = length(fx); Nx = length(x0);
if Nf ~= Nx, error('Incompatible dimensions of f and x0!'); end
if nargin < 4, MaxIter = 10; end
if nargin < 3, TolX = EPS; end
xx(1,:) = x0(:).';
%% this part is for handling Jocobina matrix
for k = 1: MaxIter

    dx = -jacobian(f,xx(k,:),h,varargin{:})\fx(:);

    xx(k + 1,:) = xx(k,:) + dx.';
    fx = feval(f,xx(k + 1,:),varargin{:});
    fxn = norm(fx);
    if fxn < TolFun || norm(dx) < TolX, break; end

end

[~,Landa_W,Landa_O] = feval(f,xx(k + 1,:),varargin{:});
x = xx(k + 1,:);
if k == MaxIter, fprintf('The best in %d iterations\n',MaxIter), end



function g = jacobian(f,x,h,varargin) 
%% this functin is for calculating the Jacobian matrix.
g=zeros(size(x,2));

if nargin < 3, h = 1e-4; end
h2 = 2*h; N = length(x); x = x(:).'; I = eye(N);

for n = 1:N
    g_plus = feval(f,x + I(n,:)*h,varargin{:});
    g_minus = feval(f,x - I(n,:)*h,varargin{:});
    g(:,n) = (g_plus - g_minus)'/h2; 

end


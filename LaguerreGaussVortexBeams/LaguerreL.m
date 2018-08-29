function v = LaguerreL(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Plot Generalized Legendre Polynomials: LaguerreL(n,z)
%
%               Coded by Manuel Diaz, NHRI, 2018.08.28.
%                   Copyright (c) 2018, Manuel Diaz.
%                           All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs.:
% [1] Milton Abramowitz, Irene Stegun, Handbook of Mathematical Functions,
%     National Bureau of Standards, 1964, ISBN: 0-486-61272-4, LC:QA47.A34.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates the Generalized Laguerre polynomials GL(n,k,x).
%
%  First terms:
%
%    GL(0,k,x) = 1;
%    GL(1,k,x) = 1 + k - x;
%    GL(2,k,x) = 2 + 3*k + k^2 - 4*x - 2*k*x + x^2)/2;
%    GL(3,k,x) = 6 + 11*k + 6*k^2 + k^3 - 18*x - 15*k*x - 3*k^2*x + ...
%                  9*x^2 + 3*k*x^2 - x^3)/6.
%  Recursion:
%
%    GL(0,k,X) = 1 
%    GL(1,k,X) = 1 + k - x;
%
%    if 2 <= a:
%
%    GL(n,k,X) = ( (k+2*n-1-X) * GL(n-1,k,X) + (1-k-n) * GL(n-2,k,X) ) / n
%
%  Special values:
%
%    For k = 0, the associated Laguerre polynomials GL(N,K,X) are equal 
%    to the Laguerre polynomials L(N,X).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program and its subprograms may be freely used, modified and
% distributed under the GNU General Public License: 
% http://www.gnu.org/copyleft/gpl.html   
%
% Basic Assumptions:
% n: is a positive interger value.
% k: is a positive interger value.
% k: can be a single value or vector array.   
%
% Function inputs:
if (nargin == 2)     % Evaluate classical Laguerre Polynomials
    n=varargin{1};
    k=0;
    x=varargin{2};
elseif (nargin == 3) % Evaluate generalized Laguerre Polynomials
    n=varargin{1};
    k=varargin{2};
    x=varargin{3};
else 
    error('Usage: >> LaguerreL(n:int,k:int,x:array)');
end

% Verify inputs
if rem(n,1)~=0, error('n must be an integer.'); end
if rem(k,1)~=0, error('k must be an integer.'); end
if n < 0, error('n must be positive integer.'); end
if k < 0, error('k must be positive integer.'); end

% Initialize solution array
v = zeros(size(x));

% Compute Laguerre Polynomials
GL = zeros( numel(x), n+1 );
if n==0
    v(:) = 1.0;         % GL(0,k,x)
elseif n==1
    v(:) = 1+k-x(:);    % GL(1,k,x)
elseif n>1
    GL(:,1) = 1;        % GL(0,k,x)
    GL(:,2) = 1+k-x(:); % GL(1,k,x)
    for i = 2:n
        GL(:,i+1) = ( (k+2*i-1-x(:)).*GL(:,i) + (-k-i+1).*GL(:,i-1) )/i;
    end
    v(:) = GL(:,i+1);   % GL(n,k,x)
end
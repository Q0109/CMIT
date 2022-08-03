function M=Matrix_Gen0(alpha,N,h,typ0)

switch typ0
    case 1   %'symToeplitz'
M=symToeplitz(alpha,N,h);
    case 2   %'LeftRL'
M=LeftRL(alpha,N,h);
    case 3   %'RightRL'
M=RightRL(alpha,N,h);
    case 5    %'RanOrt'
M=QrtCenterRiesz(alpha,N,h);
    otherwise
  error ('Matrix Generation: Wrong type of input parameters')
end
end

function R=QrtCenterRiesz(alpha,N,h)
% Make matrix R_{N}^{alpha} that corresponds
% to symmetric Riesz fractional differentiation
% using Ortiguieira's centred fractional differences
% Parameters:
%    alpha - order of differentiation (real, not necessarily integer)
%    N     - size of the resulting matrix R (N x N)
%    h     - step of discretization; default is h=1
%
% (C) 2008 Igor Podlubny, Blas Vinagre, Tomas Skovranek
%
% See:
% [1] I. Podlubny, A.Chechkin, T. Skovranek, YQ Chen,
%     B. M. Vinagre Jara, "Matrix approach to discrete
%     fractional calculus II: partial fractional differential
%     equations". http://arxiv.org/abs/0811.1355
% [2] M. Ortigueira, "Riesz potential operators and inverses
%     via fractional centred derivatives",
%     International Journal of Mathematics and Mathematical
%     Sciences, Vol. 2006, Article ID 48391, Pages 1-12,
%     DOI 10.1155/IJMMS/2006/48391

if nargin <= 1 || nargin > 3
    error('RANORT: Wrong number of input parameters')
else
    k = 0:N-1;
    rc = ((-1)*ones(size(k))).^(k).*gamma(alpha+1).*(gamma(alpha*0.5 - k + 1)).*gamma(alpha*0.5 + k + 1).^(-1);
    rc = rc * (cos(alpha*pi*0.5));
    R = zeros(N,N);
    for m = 1:N,
        R(m, m:N) = rc(1:N-m+1);
    end
    for i=1:N-1,
        for j = i:N,
            R(j,i)=R(i,j);
        end
    end
end

if nargin == 3
    R = R*h^(-alpha);
end
end

function F=RightRL(alpha,N,h)
% Make matrix F_{N}^{alpha} that corresponds
% to right-sided Riemann-Liouville fractional differentiation
% Parameters:
%    alpha - order of differentiation (real, not necessarily integer)
%    N     - size of the resulting matrix B (N x N)
%    h     - step of discretization; default is h=1
%
% (C) 2000 Igor Podlubny
%
% See the following articles:
% [1] I. Podlubny, "Matrix approach to discrete fractional calculus",
%     Fractional Calculus and Applied Analysis,
%     vol. 3, no. 4, 2000, pp. 359-386.
%     http://people.tuke.sk/igor.podlubny/pspdf/ma2dfc.pdf
% [2] I. Podlubny, A.Chechkin, T. Skovranek, YQ Chen,
%     B. M. Vinagre Jara, "Matrix approach to discrete
%     fractional calculus II: partial fractional differential
%     equations". http://arxiv.org/abs/0811.1355


n = ceil(alpha);
F = zeros(N,N);

if nargin <= 1 || nargin > 3
    error('FAN: Wrong number of input parameters')
else
    bc=fliplr(bcrecur(alpha,N-1));
    for k=1:N
       F(k,1:k)=bc((N-k+1):N);
    end
end

if nargin == 3
    F=h^(-alpha)*F;
end

F = (-1)^n*F';
end

function R=symToeplitz(alpha,N,h)
% Make matrix R_{N}^{alpha} that corresponds
% to symmetric Riesz fractional derivative as
% a half-sum of left- and right-sided Caputo derivatives.
% Parameters:
%    alpha - order of differentiation (real, not necessarily integer)
%    N     - size of the resulting matrix B (N x N)
%    h     - step of discretization; default is h=1
%
% (C) 2008 Igor Podlubny, Blas Vinagre, Tomas Skovranek
%
% See the following articles:
% [1] I. Podlubny, "Matrix approach to discrete fractional calculus",
%     Fractional Calculus and Applied Analysis,
%     vol. 3, no. 4, 2000, pp. 359-386.
%     http://people.tuke.sk/igor.podlubny/pspdf/ma2dfc.pdf
% [2] I. Podlubny, A.Chechkin, T. Skovranek, YQ Chen,
%     B. M. Vinagre Jara, "Matrix approach to discrete
%     fractional calculus II: partial fractional differential
%     equations". http://arxiv.org/abs/0811.1355

if nargin <= 1 || nargin > 3
    error ('RANSYM: Wrong number of input parameters')
else
    B = LeftRL (alpha, N+1);
    BM = B(2:(N+1), 1:N);
    R = 0.5 * (BM + BM');
end

if nargin == 3
    R = h^(-alpha)*R;
end
end

function B = LeftRL(alpha,N,h)
% Make matrix B_{N}^{alpha} that corresponds
% to left-sided Riemann-Liouville fractional differentiation
% Parameters:
%    alpha - order of differentiation (real, not necessarily integer)
%    N     - size of the resulting matrix B (N x N)
%    h     - step of discretization; default is h=1
%
% (C) 2000 Igor Podlubny
%
% See the following articles:
% [1] I. Podlubny, "Matrix approach to discrete fractional calculus",
%     Fractional Calculus and Applied Analysis,
%     vol. 3, no. 4, 2000, pp. 359-386.
%     http://people.tuke.sk/igor.podlubny/pspdf/ma2dfc.pdf
% [2] I. Podlubny, A.Chechkin, T. Skovranek, YQ Chen,
%     B. M. Vinagre Jara, "Matrix approach to discrete
%     fractional calculus II: partial fractional differential
%     equations". http://arxiv.org/abs/0811.1355

B = zeros(N,N);
if nargin <= 1 || nargin > 3
    error('BAN: Wrong number of input parameters')
else
    bc=fliplr(bcrecur(alpha,N-1));
    for k=1:N
       B(k,1:k)=bc((N-k+1):N);
    end
end

if nargin == 3
    B=h^(-alpha)*B;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=bcrecur(a, n)
%
% Computation of the fractional difference coefficients
% by the recurrence relation
%           bc(j)=(1-(a+1)/j)*bc(j-1),
%  a - order of the fractional difference
%  n - required number of coefficients
%
%  Computation of bcrecur(k,k) takes (3*k+1) flops.
%
%  Copyright (C) Igor Podlubny
%  15 Nov 1994
%
%  See also:
%  [1] Podlubny, I.: Fractional Differential Equations.
%      Academic Press, San Diego, 1999, 368 pages, ISBN 0125588402.
y=cumprod([1, 1 - ((a+1) ./ (1:n))]);
end

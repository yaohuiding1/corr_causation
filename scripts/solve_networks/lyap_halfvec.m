function out = lyap_halfvec(A, Sigma_W)
%LYAP_HALFVEC  Steady-state covariance of a linear SDE via half-vectorization.
%
%   out = LYAP_HALFVEC(A, Sigma_W) solves the continuous Lyapunov equation
%
%       A*Sigma_Y + Sigma_Y*A.' + Sigma_W = 0
%
%   for the steady-state covariance Sigma_Y of the linear stochastic system
%   dY = A Y dt + dW,  where Sigma_W is the (symmetric, typically diagonal)
%   noise covariance.
%
%   DUAL MODE
%   ---------
%   * NUMERIC input (A and Sigma_W both numeric):
%       Defers to MATLAB's built-in solver. Returns
%           out.Sigma_Y   = lyap(A, Sigma_W)        (numeric covariance)
%           out.Rho       = corr matrix of Sigma_Y  (rho_ij = s_ij/sqrt(s_ii s_jj))
%
%   * SYMBOLIC input (A and Sigma_W are sym):
%       Solves symbolically via the half-vectorization technique:
%           vec(A*S + S*A.') = (I (x) A + A (x) I) vec(S)            [Kronecker]
%           vec(S) = D_n * vech(S),   vech(S) = L_n * vec(S)         [dup/elim]
%           K = L_n (I (x) A + A (x) I) D_n
%           vech(Sigma_Y) = -K^{-1} vech(Sigma_W)
%       with the lower-triangular, column-major vech ordering
%           (1,1),(2,1),...,(n,1),(2,2),(3,2),...,(n,n),
%       and additionally records polynomial-complexity metadata (see below).
%
%   INPUTS
%       A        n-by-n causal-influence matrix (a_ij = effect of j on i).
%       Sigma_W  n-by-n noise covariance (symmetric; usually diagonal).
%
%   OUTPUT (struct) -- SYMBOLIC MODE
%       out.A            the input A (echoed for provenance)
%       out.Sigma_W      the input Sigma_W
%       out.n            system size
%       out.K            the half-vec Lyapunov operator L_n(I@A+A@I)D_n
%       out.detK         |K|  (the determinant; the common denominator)
%       out.Sigma_Y      symbolic steady-state covariance (n-by-n, simplified)
%       out.Sigma_Y_num  expand(-detK * Sigma_Y): polynomial numerators, so that
%                        each entry obeys Sigma_Y(i,j) = -Sigma_Y_num(i,j)/detK,
%                        i.e. sigma_ij = -1/|K| * Sigma_Y_num(i,j). This matches
%                        the appendices' convention exactly (verified against
%                        Appendix C/D), so Sigma_Y_num can be pasted directly
%                        as the bracketed term next to a -1/|K| prefactor with
%                        no manual sign reconciliation.
%       out.nterms       n-by-n: number of distinct monomial terms in each
%                        numerator polynomial Sigma_Y_num(i,j) (after expand).
%       out.degree       n-by-n: total polynomial degree of each numerator.
%       out.detK_nterms  number of monomial terms in |K|.
%       out.detK_degree  total polynomial degree of |K|.
%       (The nterms/degree fields summarise the combinatorial size of each
%        closed-form entry -- useful for reasoning about downstream cost.)
%
%   OUTPUT (struct) -- NUMERIC MODE
%       out.Sigma_Y      numeric steady-state covariance, from lyap().
%       out.Rho          numeric correlation matrix (Sigma_Y normalized).
%
%   See also LYAP, KRON, COEFFS, POLYNOMIALDEGREE.

    % ----- size & basic checks -----------------------------------------
    [n, nc] = size(A);
    if n ~= nc
        error('lyap_halfvec:Asquare', 'A must be square.');
    end
    if ~isequal(size(Sigma_W), [n n])
        error('lyap_halfvec:Wsize', 'Sigma_W must be the same size as A.');
    end
    if isa(A, 'sym') ~= isa(Sigma_W, 'sym')
        error('lyap_halfvec:MixedMode', ...
            'A and Sigma_W must both be symbolic or both be numeric (mixed input not supported).');
    end

    % ----- NUMERIC MODE: defer to built-in solver ----------------------
    if ~isa(A, 'sym') && ~isa(Sigma_W, 'sym')
        out.Sigma_Y = lyap(A, Sigma_W);   % solves A*X + X*A' + Sigma_W = 0
        d = sqrt(diag(out.Sigma_Y));      % std devs
        out.Rho = out.Sigma_Y ./ (d * d.');  % correlation: rho_ij = s_ij/sqrt(s_ii s_jj)
        return;
    end

    % ===================================================================
    % SYMBOLIC MODE
    % ===================================================================

    % ----- elimination / duplication matrices --------------------------
    L = elimination_matrix(n);     % m-by-n^2 , m = n(n+1)/2
    D = duplication_matrix(n);     % n^2-by-m

    % sanity check: eliminating then re-duplicating must recover vech, i.e.
    % L*D = I_m  with m = n(n+1)/2  (for n=3 this is I_6). If not, one or both
    % of L, D is malformed and every result below would be wrong.
    m = n*(n+1)/2;
    if ~isequal(L*D, eye(m))
        error('lyap_halfvec:LDcheck', ...
            'L_n*D_n ~= I_%d ; elimination/duplication matrix construction is wrong.', m);
    end

    % ----- half-vec Lyapunov operator ----------------------------------
    I = eye(n);
    T = kron(I, A) + kron(A, I);   % n^2-by-n^2
    K = L * T * D;                 % m-by-m
    K = simplify(K);

    detK = simplify(det(K));

    % ----- solve  vech(Sigma_Y) = -K^{-1} vech(Sigma_W) ----------------
    vechW = L * Sigma_W(:);                 % vech(Sigma_W) (Sigma_W symmetric)
    vechY = -(K \ vechW);                   % symbolic linear solve
    Sigma_Y = reshape(D * vechY, n, n);     % rebuild full symmetric matrix
    Sigma_Y = simplify(Sigma_Y);

    % numerator form: each entry = -Sigma_Y_num(i,j) / detK, i.e.
    % sigma_ij = -1/|K| * Sigma_Y_num(i,j) -- matches the appendices'
    % -1/|K| convention directly (verified against Appendix C/D), so no
    % manual sign reconciliation is needed when transcribing this output.
    % simplifyFraction first clears any apparent denominator -- for dense
    % topologies detK is a non-minimal common denominator (e.g. detK = 4*den),
    % so detK*Sigma_Y is a polynomial that `expand` alone may leave in rational
    % form; simplifyFraction reduces it to a genuine polynomial before expand.
    Sigma_Y_num = expand(simplifyFraction(-detK * Sigma_Y));

    % ----- polynomial-complexity metadata ------------------------------
    nterms = zeros(n);
    degree = zeros(n);
    for i = 1:n
        for j = 1:n
            [nterms(i,j), degree(i,j)] = poly_meta(Sigma_Y_num(i,j));
        end
    end
    [detK_nterms, detK_degree] = poly_meta(detK);

    % ----- pack output --------------------------------------------------
    out.A           = A;
    out.Sigma_W     = Sigma_W;
    out.n           = n;
    out.K           = K;
    out.detK        = detK;
    out.Sigma_Y     = Sigma_Y;
    out.Sigma_Y_num = Sigma_Y_num;
    out.nterms      = nterms;
    out.degree      = degree;
    out.detK_nterms = detK_nterms;
    out.detK_degree = detK_degree;
end


% ======================================================================
% local helper functions
% ======================================================================

function L = elimination_matrix(n)
%ELIMINATION_MATRIX  L_n : vech(S) = L_n vec(S), lower-tri column-major.
    m = n*(n+1)/2;
    L = zeros(m, n^2);
    k = 0;
    for j = 1:n          % column
        for i = j:n      % row (lower triangle)
            k = k + 1;
            col = (j-1)*n + i;     % column-major vec index of (i,j)
            L(k, col) = 1;
        end
    end
end


function D = duplication_matrix(n)
%DUPLICATION_MATRIX  D_n : vec(S) = D_n vech(S) for symmetric S.
    m = n*(n+1)/2;
    D = zeros(n^2, m);
    for j = 1:n
        for i = 1:n
            row = (j-1)*n + i;             % vec index of (i,j)
            ii  = max(i,j);  jj = min(i,j);
            k   = vech_index(ii, jj, n);   % vech index of (max,min)
            D(row, k) = 1;
        end
    end
end


function k = vech_index(i, j, n)
%VECH_INDEX  position of lower-tri entry (i,j), i>=j, in column-major vech.
    k = 0;
    for c = 1:j-1
        k = k + (n - c + 1);   % entries in columns before column j
    end
    k = k + (i - j + 1);       % offset within column j
end


function [nt, dg] = poly_meta(p)
%POLY_META  number of monomial terms and total degree of a polynomial.
    p = numden(p);          % robust: take the numerator if p is a rational expr
    p = expand(p);
    if isequal(p, sym(0))
        nt = 0;  dg = 0;  return;     % zero polynomial
    end
    v = symvar(p);
    if isempty(v)
        nt = 1;  dg = 0;  return;     % nonzero constant
    end
    c  = coeffs(p, v);                 % one coefficient per distinct monomial
    nt = numel(c);
    dg = double(polynomialDegree(p));  % total degree across all variables
end

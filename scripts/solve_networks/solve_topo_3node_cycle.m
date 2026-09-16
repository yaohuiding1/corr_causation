%% solve_topo_3node_cycle.m
% Topology: 3-node directed cycle, 1->3->2->1 (edges a31, a23, a12).
%
% Convention: a_ij = causal influence of node j on node i (j -> i).
% Diagonal generalized to distinct a11, a22, a33.
% NOTE: the special cyclic parameterization (e.g. a12 = -a23*a31) is NOT
% imposed here; that constrained case is left for separate study.
%
% Solves  A*Sigma_Y + Sigma_Y*A' + Sigma_W = 0  via lyap_halfvec.m.

clc; clear;

% --- symbols ---------------------------------------------------------
syms a11 a22 a33 a12 a23 a31
syms sigmaW1 sigmaW2 sigmaW3

% --- causal matrix A  (a_ij = effect of j on i) ----------------------
% Off-diagonals: a12 (2->1), a23 (3->2), a31 (1->3).
A = [a11  a12    0 ;
       0  a22  a23;
     a31    0  a33];

% --- noise covariance (diagonal, independent noise) ------------------
Sigma_W = diag([sigmaW1 sigmaW2 sigmaW3]);

% --- node-permutation symmetry of the wiring -------------------------
network_sym = network_symmetry(A);
disp(['Network symmetry: ' network_sym.label]);

% --- solve -----------------------------------------------------------
out = lyap_halfvec(A, Sigma_W);

% --- report ----------------------------------------------------------
disp('|K| (determinant) =');              disp(out.detK);
disp('Sigma_Y numerators (entry/|K|) =');  disp(out.Sigma_Y_num);
disp('# terms per numerator =');           disp(out.nterms);
disp('total degree per numerator =');      disp(out.degree);

% headline correlation for this topology: rho_21
rho21 = simplify(out.Sigma_Y(2,1) / sqrt(out.Sigma_Y(1,1) * out.Sigma_Y(2,2)));
disp('rho_21 ='); disp(rho21);

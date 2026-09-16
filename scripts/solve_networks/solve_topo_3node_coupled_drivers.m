%% solve_topo_3node_coupled_drivers.m
% Topology: 3-node coupled drivers. Nodes 1 and 2 are reciprocally coupled
% (1<->2) and both drive node 3 (1->3, 2->3).
%
% Convention: a_ij = causal influence of node j on node i (j -> i).
% The diagonal is generalized to distinct a11, a22, a33 here so that
% different parameterizations (e.g. setting them equal) can be tried below.
%
% Solves  A*Sigma_Y + Sigma_Y*A' + Sigma_W = 0  via lyap_halfvec.m.

clc; clear;

% --- symbols ---------------------------------------------------------
syms a11 a22 a33 a12 a21 a31 a32
syms sigmaW1 sigmaW2 sigmaW3

% --- causal matrix A  (a_ij = effect of j on i) ----------------------
% Off-diagonals: a12 (2->1), a21 (1->2), a31 (1->3), a32 (2->3).
% Edit this block to experiment with other parameterizations.
A = [a11   a12    0 ;
     a21   a22    0 ;
     a31   a32   a33];

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

% headline correlation for this topology: rho_21 (between the two drivers)
rho21 = simplify(out.Sigma_Y(2,1) / sqrt(out.Sigma_Y(1,1) * out.Sigma_Y(2,2)));
disp('rho_21 ='); disp(rho21);

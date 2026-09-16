%% solve_topo_2node_dyad.m
% Topology: 2-node dyadic (reciprocal) interaction, 1<->2.
%
% Convention: a_ij = causal influence of node j on node i (j -> i).
% Diagonal generalized to distinct a11, a22 (edit A below to experiment).
%
% Solves  A*Sigma_Y + Sigma_Y*A' + Sigma_W = 0  via lyap_halfvec.m.

clc; clear;

% --- symbols ---------------------------------------------------------
syms a11 a22 a12 a21
syms sigmaW1 sigmaW2

% --- causal matrix A  (a_ij = effect of j on i) ----------------------
% Off-diagonals: a12 (2->1), a21 (1->2).
A = [a11  a12;
     a21  a22];

% --- noise covariance (diagonal, independent noise) ------------------
Sigma_W = diag([sigmaW1 sigmaW2]);

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

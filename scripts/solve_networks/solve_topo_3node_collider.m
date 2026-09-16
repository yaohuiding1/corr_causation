%% solve_topo_3node_collider.m
% Topology: 3-node collider. Nodes 1 and 2 both drive node 3 (1->3, 2->3);
% no link between 1 and 2.
%
% Convention: a_ij = causal influence of node j on node i (j -> i).
% Diagonal generalized to distinct a11, a22, a33.
%
% Solves  A*Sigma_Y + Sigma_Y*A' + Sigma_W = 0  via lyap_halfvec.m.

clc; clear;

% --- symbols ---------------------------------------------------------
syms a11 a22 a33 a31 a32
syms sigmaW1 sigmaW2 sigmaW3

% --- causal matrix A  (a_ij = effect of j on i) ----------------------
% Off-diagonals: a31 (1->3), a32 (2->3).
A = [a11   0    0 ;
       0  a22   0 ;
     a31  a32  a33];

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

% headline correlations: rho_21 (the two causes; = 0 for a pure collider)
% and rho_31 (a cause and the collider).
rho21 = simplify(out.Sigma_Y(2,1) / sqrt(out.Sigma_Y(1,1) * out.Sigma_Y(2,2)));
rho31 = simplify(out.Sigma_Y(3,1) / sqrt(out.Sigma_Y(1,1) * out.Sigma_Y(3,3)));
disp('rho_21 ='); disp(rho21);
disp('rho_31 ='); disp(rho31);

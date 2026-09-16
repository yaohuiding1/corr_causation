function netsym = network_symmetry(A)
%NETWORK_SYMMETRY  Node-permutation (automorphism) symmetry of a network.
%
%   netsym = NETWORK_SYMMETRY(A) returns the group of node relabelings that
%   leave the wiring of A unchanged. A permutation p of the nodes is a
%   symmetry (graph automorphism) iff it preserves the sparsity pattern of
%   A, i.e.  A(i,j) ~= 0  <=>  A(p(i),p(j)) ~= 0  for all i,j.
%
%   Diagonal/self-loops are always present and map among themselves, so they
%   never affect the result; only the off-diagonal pattern matters.
%
%   For n nodes all n! permutations are tested (cheap for n <= 4). The full
%   group captures BOTH reflections (transpositions, e.g. "labels 1 and 2 are
%   interchangeable") AND rotations (cyclic relabelings) -- for 3 nodes this is
%   the dihedral symmetry of a triangle, S_3 = D_3.
%
%   INPUT
%       A    n-by-n causal matrix (symbolic or numeric).
%
%   OUTPUT (struct)
%       netsym.perms     k-by-n matrix; each row is a symmetry permutation p
%                     (p(i) = new label of node i), identity always included.
%       netsym.order     k = number of symmetries (1 means none/trivial).
%       netsym.label     human-readable description.
%       netsym.hasSymmetry  true if order > 1.
%
%   Examples (3-node):
%       fully connected -> full S_3 (all 6 permutations)
%       coupled drivers -> reflection (1 2)
%       directed cycle  -> rotational (1 2 3),(1 3 2)
%       mediator        -> trivial (no symmetry)
%
%   See also PERMS.

    n = size(A,1);

    % ----- sparsity pattern (logical) ----------------------------------
    if isa(A, 'sym')
        P = ~isAlways(A == 0, 'Unknown', 'false');   % true where entry is nonzero
    else
        P = (A ~= 0);
    end

    % ----- brute-force all n! permutations -----------------------------
    allP = perms(1:n);                  % each row a permutation
    keep = false(size(allP,1),1);
    for r = 1:size(allP,1)
        p = allP(r,:);
        if isequal(P(p,p), P)           % pattern preserved under relabeling
            keep(r) = true;
        end
    end
    G = allP(keep,:);

    % put identity first, then sort the rest for stable output
    idRow = all(G == (1:n), 2);
    G = [G(idRow,:); sortrows(G(~idRow,:))];

    % ----- readable label ----------------------------------------------
    k = size(G,1);
    netsym.perms      = G;
    netsym.order      = k;
    netsym.hasSymmetry = k > 1;

    if k == 1
        netsym.label = 'trivial -- no node symmetry';
        return;
    end

    % describe non-identity elements in cycle notation
    cyc = strings(0,1);
    for r = 1:k
        if ~all(G(r,:) == (1:n))
            cyc(end+1,1) = string(perm_cycles(G(r,:))); %#ok<AGROW>
        end
    end
    cycStr = char(strjoin(cyc, ", "));

    isTransp = @(p) nnz(p ~= (1:n)) == 2;   % swaps exactly two nodes
    nonId = G(~all(G == (1:n),2), :);
    allTransp = all(arrayfun(@(r) isTransp(nonId(r,:)), 1:size(nonId,1)));

    if k == factorial(n)
        netsym.label = sprintf('full S_%d (all node labels interchangeable): %s', n, cycStr);
    elseif k == 2 && allTransp
        netsym.label = sprintf('reflection %s', cycStr);
    elseif allTransp
        netsym.label = sprintf('reflections %s', cycStr);
    elseif ~any(arrayfun(@(r) isTransp(nonId(r,:)), 1:size(nonId,1)))
        netsym.label = sprintf('rotational (cyclic) %s', cycStr);
    else
        netsym.label = sprintf('order-%d group: %s', k, cycStr);
    end
end


function s = perm_cycles(p)
%PERM_CYCLES  cycle notation of a permutation (1-line, identity -> "e").
    n = numel(p);
    seen = false(1,n);
    parts = strings(0,1);
    for i = 1:n
        if ~seen(i)
            c = i; seen(i) = true; nxt = p(i);
            while nxt ~= i
                c(end+1) = nxt; seen(nxt) = true; nxt = p(nxt); %#ok<AGROW>
            end
            if numel(c) > 1
                parts(end+1,1) = "(" + strjoin(string(c), " ") + ")"; %#ok<AGROW>
            end
        end
    end
    if isempty(parts), s = "e"; else, s = strjoin(parts, ""); end
end

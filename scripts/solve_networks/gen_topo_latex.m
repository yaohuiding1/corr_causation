% gen_topo_latex.m
% Emits ready-to-paste LaTeX (Determinant + covariance, "-1/|K| factored once",
% grouped by noise source) for each network topology, straight from the
% symbolic solution produced by lyap_halfvec. NOT inserted automatically --
% the printed blocks are copied by hand into appendices/appendix_d_triadic_interaction.tex
%
% Each covariance entry is split as  num = P1*sigmaW1 + P2*sigmaW2 + P3*sigmaW3
% (covariance is exactly linear in the noise variances), and rendered as
%   sigma_ij = (-1/|K|)( P1 sigma_{W1} + P2 sigma_{W2} + P3 sigma_{W3} ).
clc; clear;

names = {'dyad','confounder','mediator','collider','moderator', ...
         'coupled_drivers','cycle','mediator_2way','full'};

for t = 1:numel(names)
    name = names{t};
    fprintf('\n\n======================== %s ========================\n', name);
    emit_topo(name);
end


% ----------------------------------------------------------------------
function emit_topo(name)
    syms a11 a22 a33 a12 a21 a13 a31 a23 a32 a real
    syms sigmaW1 sigmaW2 sigmaW3 real

    [A, Sw, subEqual] = topo_def(name);
    n = size(A,1);

    out = lyap_halfvec(A, Sw);
    detK = out.detK;
    num  = out.Sigma_Y_num;

    if subEqual
        detK = simplify(subs(detK, [a11 a22 a33], [a a a]));
        num  = simplify(subs(num,  [a11 a22 a33], [a a a]));
    end

    % ---- determinant, expanded (not factored -- |K| is mostly bookkeeping,
    % rarely makes it into the main manuscript, so a factored form isn't worth
    % the risk of MATLAB's factor() landing on an awkward/meaningless split) ----
    fprintf('|\\mathbf{K}| = %s\n', latex(expand(detK)));

    % ---- covariance entries, lower triangular column-major ----
    nz = [sigmaW1 sigmaW2 sigmaW3];
    nz = nz(1:n);
    for j = 1:n
        for i = j:n
            e = expand(num(i,j));
            % split by noise source -> one line per nonzero noise group
            pieces = {};
            for k = 1:n
                unit = zeros(1,n); unit(k) = 1;
                Pk = expand(subs(e, nz, unit));
                if isequal(Pk, sym(0)); continue; end
                coef = latex(Pk);
                % parenthesize multi-term coefficients so sigma_{Wk} binds to all
                bare = coef; if startsWith(bare,'-'); bare = bare(2:end); end
                if contains(bare,'+') || contains(bare,'-')
                    coef = ['\left(' coef '\right)'];
                end
                term = sprintf('%s\\,\\sigma_{W_{%d}}', coef, k);
                if isempty(pieces)
                    pieces{end+1} = term;                 % first group
                elseif startsWith(coef,'-')
                    pieces{end+1} = ['{}' term];          % continuation: binary minus
                else
                    pieces{end+1} = ['{}+ ' term];        % continuation: binary plus
                end
            end
            % consistency check: numerator == sum_k Pk sigmaWk
            chk = e;
            for k = 1:n
                unit = zeros(1,n); unit(k) = 1;
                chk = chk - expand(subs(e, nz, unit))*nz(k);
            end
            if ~isequal(expand(chk), sym(0))
                fprintf('%% WARNING: sigma_%d%d not linear in noise!\n', i, j);
            end
            if isempty(pieces); pieces = {'0'}; end
            if numel(pieces) == 1
                fprintf('\\sigma_{%d%d}&= -\\frac{1}{|\\mathbf{K}|}\\left( %s \\right)\\\\\n', ...
                        i, j, pieces{1});
            else
                fprintf('\\sigma_{%d%d}&= -\\frac{1}{|\\mathbf{K}|}\\big( %s\\\\\n', ...
                        i, j, pieces{1});
                for p = 2:numel(pieces)-1
                    fprintf('&\\qquad %s\\\\\n', pieces{p});
                end
                fprintf('&\\qquad %s \\big)\\\\\n', pieces{end});
            end
        end
    end
end


% ----------------------------------------------------------------------
function [A, Sw, subEqual] = topo_def(name)
    syms a11 a22 a33 a12 a21 a13 a31 a23 a32 real
    syms sigmaW1 sigmaW2 sigmaW3 real
    switch name
        case 'dyad'                       % 2-node, distinct diagonal
            A = [a11 a12; a21 a22];
            Sw = diag([sigmaW1 sigmaW2]);   subEqual = false;
        case 'confounder'                   % distinct diagonal
            A = [a11 0 a13; 0 a22 a23; 0 0 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = false;
        case 'mediator'                     % single a
            A = [a11 0 0; 0 a22 a23; a31 0 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = true;
        case 'collider'                     % single a; 1->3, 2->3
            A = [a11 0 0; 0 a22 0; a31 a32 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = true;
        case 'moderator'                    % single a
            A = [a11 a12 a13; a21 a22 a23; 0 0 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = true;
        case 'coupled_drivers'              % single a
            A = [a11 a12 0; a21 a22 0; a31 a32 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = true;
        case 'cycle'                        % single a
            A = [a11 a12 0; 0 a22 a23; a31 0 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = true;
        case 'mediator_2way'               % single a
            A = [a11 0 a13; 0 a22 a23; a31 a32 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = true;
        case 'full'                         % distinct diagonal
            A = [a11 a12 a13; a21 a22 a23; a31 a32 a33];
            Sw = diag([sigmaW1 sigmaW2 sigmaW3]); subEqual = false;
        otherwise
            error('unknown topology %s', name);
    end
end

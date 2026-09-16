%% Contour + vector-field plot of rho_21 for the 3-node moderator topology,
% as a function of a12, a21 (see figure caption for what each plotted
% element represents). Closed-form rho_21 derived in Appendix D
% (eq:rho21_moderator / main text eq:rho21_moderator_main).
%
% One standalone square PDF per panel, plus one standalone colorbar PDF,
% assembled in LaTeX per the authoritative guide
% (tasks/stage3_tasks_plot_networks/stage3_topo_x-node_figure_template_instructions.md).
% Structure copied from plot_topo_3node_cycle_subplot.m (the reference for
% topologies WITH a stability mask), since the moderator's dyadic-block
% condition (a12*a21 < a^2) is a genuine mask, not the acyclic no-mask case.
% Output: figures/topo_3node_moderator_figures/{B1,B2,B3,C1,C2,C3,colorbar}.pdf
%
% PANEL MAPPING: a, sigma_W1, sigma_W3 are fixed across all six panels;
% (a13,a23) varies per panel, and sigma_W2 varies by row (row B = sigma_W2=1,
% row C = sigma_W2=3), giving a row/column contrast similar to
% confounder/collider/mediator's sigma_W3 row split.
clc; clear;

a = -0.5;
sigma_W1 = 1; sigma_W3 = 1;
% sigma_W2 differs by row: row B (panels 1-3) keeps sigma_W2=1, row C
% (panels 4-6) uses sigma_W2=3, to explore variance asymmetry between
% nodes 1 and 2.
sigma_W2_by_row = [1, 1, 1, 3, 3, 3];

% Six (a13,a23) combinations, none the antipode of another: since rho_21
% depends on (a13,a23) only through a13^2, a23^2, and a13*a23 (node 3 has
% no feedback from nodes 1/2), sign-negated pairs (e.g. (-0.5,-0.5) and
% (0.5,0.5)) are mathematically forced to produce identical contours -- so
% the six values below span zero/single-driver/orthogonal/same-sign/
% asymmetric-opposite-sign cases within [-1,1] without any such redundant
% pair. Panel 1 is exactly (0,0), recovering the dyadic formula exactly
% (Appendix D).
params = [...
     0.0,    0.0;
     1.0,    0.0;
     0.0,    1.0;
     1.0,    1.0;
     0.5,   -1.0;
    -1.0,    0.5];
% Row A is reserved for the network diagram and matrices strip (drawn
% directly in LaTeX, not exported here), so the two MATLAB-generated data
% rows are B and C.
panel_labels = {'B1','B2','B3','C1','C2','C3'};

% param grid for a12 and a21. Range [-1.5,1.5], matching the stage-wide
% convention.
range = 1.5;
N = 150;
[a12, a21] = meshgrid(linspace(-range, range, N), linspace(-range, range, N));
a12sq = a12.^2;
a21sq = a21.^2;
a12a21 = a12.*a21;
a12cube = a12.^3;
a21cube = a21.^3;

%% Diverging blue-white-red colormap (stage-wide convention).
rdbu_anchors = [ ...
      5,  48,  97;   33, 102, 172;   67, 147, 195;  146, 197, 222; ...
    209, 229, 240;  247, 247, 247;  253, 219, 199;  244, 165, 130; ...
    214,  96,  77;  178,  24,  43;  103,   0,  31] / 255;
neutral_halfwidth = 0.05;
v = linspace(-1, 1, 256);
t = 0.5 + 0.5 * sign(v) .* max(abs(v) - neutral_halfwidth, 0) ...
                        ./ (1 - neutral_halfwidth);
cmap_div = interp1(linspace(0,1,size(rdbu_anchors,1)), rdbu_anchors, t);

invalid_grey = [0.80 0.80 0.80];
rho_fill_levels = -1:0.02:1;
rho_levels = linspace(-1,1,31);  % 30 visible lines (excluding zero), matching
                                  % the confounder/collider/mediator/cycle
                                  % convention

%% Stability mask: a and (a12,a21) alone determine stability; a13, a23 have
% NO effect (block-triangular A: node 3 receives no feedback). Stable iff
% a<0 (fixed here) AND a12*a21 < a^2. Confirmed via Appendix D:
% |K| = 16a^2(a^2-a12*a21)(4a^2-a12*a21) is strictly positive under Hurwitz
% stability, so no sign(K) factor is needed in rho_21 below.

%% Output folder -- one subfolder per topology.
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fullfile(script_dir, '..', '..');
fig_dir    = fullfile(proj_root, 'figures', 'topo_3node_moderator_figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Fixed square canvas and fixed axes Position -- IDENTICAL for every panel.
panel_size = 360;
ax_pos = [0.16 0.12 0.70 0.70];

%% Per-panel: compute rho_* gradient, plot, export.
for k = 1:size(params,1)
    a13 = params(k,1);
    a23 = params(k,2);
    sigma_W2 = sigma_W2_by_row(k);

    % The leading minus sign is correct: numerator/denominators carry a
    % common 4a scale factor relative to Appendix D's N21/D11/D22;
    % sqrt(denom1*denom2) picks up |4a|=-4a (a<0), so the surviving sign is
    % -sign(a)=+1, recovering Appendix D's rho_21 exactly.
    numerator = (-16*a^4*a21+4*a^2*a12.*a21sq)*sigma_W1 ...
        + (-16*a^4*a12+4*a^2*a21.*a12sq)*sigma_W2 ...
        + (-12*(a^2)*(a23^2)*a12 -12*(a^2)*(a13^2)*a21 ...
           + 8*a*a13*a23*a12a21 +16*(a^3)*a13*a23)*sigma_W3;

    denom1 = (32*(a^5)-24*(a^3)*a12a21+4*a*a12sq.*a21sq)*sigma_W1 ...
        + (16*(a^3)*a12sq - 4*a*a12cube.*a21)*sigma_W2 ...
        + (16*(a^3)*a13^2 - 24*a^2*a12*a13*a23+12*a*a12sq*a23^2 ...
           - 4*a*a12a21*a13^2)*sigma_W3;

    denom2 = (16*(a^3)*a21sq - 4*a*a21cube.*a12)*sigma_W1 ...
        + (32*(a^5)-24*(a^3)*a12a21+4*a*a12sq.*a21sq)*sigma_W2 ...
        + (16*(a^3)*a23^2 - 24*a^2*a21*a23*a13+12*a*a21sq*a13^2 ...
           - 4*a*a12a21*a23^2)*sigma_W3;

    rho_21 = -numerator ./ sqrt(denom1 .* denom2);

    % Stability mask: a<0 (fixed, true here); binding condition a12*a21<a^2.
    valid = a12a21 < a^2;
    rho_21_masked = rho_21;
    rho_21_masked(~valid) = NaN;

    [drho_da12, drho_da21] = gradient(rho_21_masked, a12(1,:), a21(:,1));
    drho_da12(~valid) = NaN;
    drho_da21(~valid) = NaN;

    fig = figure('Color','w','Units','points', ...
        'Position',[100 100 panel_size panel_size],'Visible','off');
    ax = axes(fig,'Units','normalized','Position',ax_pos);

    contourf(ax, a12, a21, rho_21_masked, rho_fill_levels, 'LineStyle', 'none')
    hold(ax,'on')
    contour(ax, a12, a21, rho_21_masked, rho_levels(rho_levels ~= 0), ...
            'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.5)
    colormap(ax, cmap_div)
    set(ax, 'Color', invalid_grey)
    clim(ax, [-1 1])

    % rho_21 = 0 contour, bright green dashed (stage-wide convention).
    contour(ax, a12, a21, rho_21_masked, [0 0], '--', ...
            'LineColor', [0.00 0.62 0.05], 'LineWidth', 2.5)
    % Stability boundary, black solid.
    contour(ax, a12, a21, a12a21, [a*a a*a], 'k-', 'LineWidth', 2.5)

    mag = sqrt(drho_da12.^2 + drho_da21.^2);
    U = drho_da12 ./ mag;
    V = drho_da21 ./ mag;
    step = 20;
    arrow_scale = 0.5;
    quiver(ax, a12(1:step:end,1:step:end), ...
        a21(1:step:end,1:step:end), ...
        U(1:step:end,1:step:end), ...
        V(1:step:end,1:step:end), ...
        arrow_scale,'k', 'LineWidth', 1)

    axis(ax,'square')

    % set FontSize before xlabel/ylabel/title (ordering pitfall).
    set(ax,'FontSize',16)

    axis_ticks = -range:0.5:range;
    set(ax, 'XTick', axis_ticks, 'YTick', axis_ticks)
    ax.XAxis.TickLabelRotation = 0;
    ax.YAxis.TickLabelRotation = 0;

    % Font sizes match the dyad/confounder/collider/mediator/two-way-mediator/
    % cycle reference: axis labels 22, title 18, ticks 16.
    xlabel(ax,'\textbf{$\mathbf{a_{12}}$}','Interpreter','latex','FontSize',22)
    ylabel(ax,'\textbf{$\mathbf{a_{21}}$}','Interpreter','latex','FontSize',22)

    % %+.2f (explicit sign) for a13/a23, matching the stage-wide sign
    % convention for parameters that vary in sign across panels.
    title(ax, sprintf('$a_{13}=%+.2f$, $a_{23}=%+.2f$', a13, a23), ...
        'Interpreter','latex','FontSize',18);

    hold(ax,'off')

    out_path = fullfile(fig_dir, [panel_labels{k} '.pdf']);
    exportgraphics(fig, out_path, 'ContentType', 'vector');
    fprintf('Wrote %s\n', out_path);
    close(fig);
end

%% Standalone colorbar (identical construction to the dyad/confounder/
% collider/mediator/two-way-mediator/cycle reference). Label is rho_21.
cb_width  = 3 * panel_size * 1.1;
cb_height = 160;
fig = figure('Color','w','Units','points','Position',[100 100 cb_width cb_height], ...
    'Visible','off');
ax = axes(fig,'Visible','off');
colormap(ax, cmap_div)
clim(ax, [-1 1])
cb = colorbar(ax,'Location','south');
cb.AxisLocation = 'out';
cb.Position = [0.08 0.35 0.84 0.20];
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathbf{\rho_{21}}$';
cb.Label.FontWeight = 'bold';
cb.FontSize = 18;
cb.Label.FontSize = 22;

cb_path = fullfile(fig_dir, 'colorbar.pdf');
exportgraphics(fig, cb_path, 'ContentType', 'vector');
fprintf('Wrote %s\n', cb_path);
close(fig);

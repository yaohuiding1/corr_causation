%% Contour + vector-field plot of rho_21 for the 3-node cycle topology, as a
% function of a31, a23 (see figure caption for what each plotted element
% represents). Closed-form rho_21 derived in Appendix D (eq:rho_cycle /
% main text eq:rho21_cycle_main).
%
% One standalone square PDF per panel, plus one standalone colorbar PDF,
% assembled in LaTeX per the authoritative guide
% (tasks/stage3_tasks_plot_networks/stage3_topo_x-node_figure_template_instructions.md).
% Structure copied from plot_topo_3node_mediator_2way_subplot.m (the
% reference for topologies WITH a stability mask), since the cycle genuinely
% needs one -- a two-sided bound, unlike the two-way mediator's one-sided one.
% Output: figures/topo_3node_cycle_figures/{B1,B2,B3,C1,C2,C3,colorbar}.pdf
%
% PANEL MAPPING: the six panels follow a per-row narrative (row 1 =
% mediator-limit baseline with a sigma_W3 contrast; row 2 = positive-quadrant
% small/moderate |a12| contrast; row 3 = negative-quadrant small/moderate
% |a12| contrast) rather than a clean 3-column single-axis contrast, so
% (like the two-way mediator) this does not transpose cleanly to the
% template's 2x3 grid. Same six parameter sets and reading order as the
% original combined figure, relabeled B1,B2,B3,C1,C2,C3 in sequence.
clc; clear;

a = -0.5;
sigma_W1 = 1;

% columns: [a12, sigma_W2, sigma_W3].
params = [...
     0.00,  1,  1;
     0.00,  3,  1;
     0.10,  1,  1;
     0.30,  1,  1;
    -0.05,  3,  1;
    -0.30,  3,  1];
% Row A is reserved for the network diagram and matrices strip (drawn
% directly in LaTeX, not exported here), so the two MATLAB-generated data
% rows are B and C.
panel_labels = {'B1','B2','B3','C1','C2','C3'};

% param grid for a31 and a23. Range [-1.5,1.5], matching the stage-wide
% convention.
range = 1.5;
N = 150;
[a31, a23] = meshgrid(linspace(-range, range, N), linspace(-range, range, N));
a31sq = a31.^2;
a23sq = a23.^2;
a23cube = a23.^3;

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
                                  % the combined script's tuned level set

%% Stability mask: TWO-SIDED bound (see
% tasks/stage3_tasks_plot_networks/stage3_plot_3node_cycle.md for the
% derivation). The cycle's characteristic polynomial is (lambda-a)^3=c with
% c=a12*a23*a31; the real root needs c<-a^3, the complex pair needs c>8a^3,
% so Hurwitz stability requires 8a^3 < c < -a^3 (a<0). On this region
% |K|=8(8a^3-c)(a^3+c)>0 identically, so no sign prefactor is applied to
% rho_21 below.
stab_lo = 8*a^3;
stab_hi = -a^3;

%% Output folder -- one subfolder per topology.
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fullfile(script_dir, '..', '..');
fig_dir    = fullfile(proj_root, 'figures', 'topo_3node_cycle_figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Fixed square canvas and fixed axes Position -- IDENTICAL for every panel.
panel_size = 360;
ax_pos = [0.16 0.12 0.70 0.70];

%% Per-panel: compute rho_* gradient, plot, export.
for k = 1:size(params,1)
    a12      = params(k,1);
    sigma_W2 = params(k,2);
    sigma_W3 = params(k,3);

    % rho_21 for the cycle (Appendix D eq:rho_cycle / main text
    % eq:rho21_cycle_main; sgn(|K|)=+1 on the stable region, see mask note).
    N21 = (-2*a^3*(a23.*a31) + a12*(a23sq.*a31sq))*sigma_W1 ...
        + (4*a^4*a12 + a*a12^2*(a23.*a31))*sigma_W2 ...
        + 3*a^2*a12*a23sq*sigma_W3;

    D11 = (8*a^5*sigma_W1 + 5*a^2*a12*sigma_W1*(a23.*a31)) ...
        + (4*a^3*a12^2 + a12^3*(a23.*a31))*sigma_W2 ...
        + 3*a*a12^2*a23sq*sigma_W3;

    D22 = (8*a^5*sigma_W2 + 5*a^2*a12*sigma_W2*(a23.*a31)) ...
        + (4*a^3*a23sq + a12*(a23cube.*a31))*sigma_W3 ...
        + 3*a*(a23sq.*a31sq)*sigma_W1;

    rho_21 = N21 ./ sqrt(D11 .* D22);

    % Stability mask: c = a12*a23*a31 must satisfy stab_lo < c < stab_hi.
    c = a12 .* a23 .* a31;
    valid = (c > stab_lo) & (c < stab_hi);
    rho_21_masked = rho_21;
    rho_21_masked(~valid) = NaN;

    [drho_da31, drho_da23] = gradient(rho_21_masked, a31(1,:), a23(:,1));
    drho_da31(~valid) = NaN;
    drho_da23(~valid) = NaN;

    fig = figure('Color','w','Units','points', ...
        'Position',[100 100 panel_size panel_size],'Visible','off');
    ax = axes(fig,'Units','normalized','Position',ax_pos);

    contourf(ax, a31, a23, rho_21_masked, rho_fill_levels, 'LineStyle', 'none')
    hold(ax,'on')
    contour(ax, a31, a23, rho_21_masked, rho_levels(rho_levels ~= 0), ...
            'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.5)
    colormap(ax, cmap_div)
    set(ax, 'Color', invalid_grey)
    clim(ax, [-1 1])

    % rho_21 = 0 contour, bright green dashed (stage-wide convention).
    contour(ax, a31, a23, rho_21_masked, [0 0], '--', ...
            'LineColor', [0.00 0.62 0.05], 'LineWidth', 2.5)

    % Stability boundary: black solid, drawn as the two level curves of the
    % field c = a12*a23*a31 at c = stab_lo and c = stab_hi. For a12=0, c is
    % identically 0 (whole plane stable) -- contour emits a harmless
    % "constant ZData" warning and draws nothing, which is correct.
    contour(ax, a31, a23, a12.*a23.*a31, [stab_lo stab_hi], 'k-', 'LineWidth', 2.5)

    mag = sqrt(drho_da31.^2 + drho_da23.^2);
    U = drho_da31 ./ mag;
    V = drho_da23 ./ mag;
    step = 20;
    arrow_scale = 0.5;
    quiver(ax, a31(1:step:end,1:step:end), ...
        a23(1:step:end,1:step:end), ...
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

    % Font sizes match the dyad/confounder/collider/mediator/two-way-mediator
    % reference: axis labels 22, title 18, ticks 16.
    xlabel(ax,'\textbf{$\mathbf{a_{31}}$}','Interpreter','latex','FontSize',22)
    ylabel(ax,'\textbf{$\mathbf{a_{23}}$}','Interpreter','latex','FontSize',22)

    % Title shows all three varying quantities (a12, sigma_W2, sigma_W3);
    % %+.2f for a12 since it varies in sign across panels, matching the
    % stage-wide sign convention.
    title(ax, sprintf('$a_{12}=%+.2f$, $\\sigma_{W_2}=%.2f$, $\\sigma_{W_3}=%.2f$', ...
        a12, sigma_W2, sigma_W3), 'Interpreter','latex','FontSize',18);

    hold(ax,'off')

    out_path = fullfile(fig_dir, [panel_labels{k} '.pdf']);
    exportgraphics(fig, out_path, 'ContentType', 'vector');
    fprintf('Wrote %s\n', out_path);
    close(fig);
end

%% Standalone colorbar (identical construction to the dyad/confounder/
% collider/mediator/two-way-mediator reference). Label is rho_21.
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

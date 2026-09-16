%% Contour + vector-field plot of rho_21 for the 3-node confounder topology,
% as a function of a13, a23 (see figure caption for what each plotted
% element represents). Closed-form rho_21 derived in Appendix D
% (eq:d26_rho21_confounder / main text eq:rho21_confounder).
%
% One standalone square PDF per panel, plus one standalone colorbar PDF,
% assembled in LaTeX per the authoritative guide
% (tasks/stage3_tasks_plot_networks/stage3_topo_x-node_figure_template_instructions.md).
% Structure copied from plot_topo_2node_dyad_subplot.m (the reference
% implementation); only the topology-specific pieces differ (formula, axes,
% the absence of a stability mask, and the param layout).
% Output: figures/topo_3node_confounder_figures/{B1,B2,B3,C1,C2,C3,colorbar}.pdf
clc; clear;

% Params for the confounder, panels B1-C3 (2x3 grid). columns: [a, sigma_W2,
% sigma_W3]; sigma_W1 = 1 throughout. Row B fixes sigma_W2=sigma_W3=1 and
% varies a in {-0.5,-1.0,-1.5} across columns; row C varies sigma_W2 and/or
% sigma_W3 per column too (see matrix values directly).
params = [...
    -0.5, 1, 1;
    -1.0, 1, 1;
    -1.5, 1, 1;
    -0.5, 3, 1;
    -1.0, 1, 3;
    -1.5, 3, 3];
sigma_W1 = 1;
% Row A is reserved for the network diagram and matrices strip (drawn
% directly in LaTeX, not exported here), so the two MATLAB-generated data
% rows are B and C.
panel_labels = {'B1','B2','B3','C1','C2','C3'};

% param grid for a13 and a23 (the two connections from node 3 into nodes 1 and 2)
range = 1.5;
N = 150;
[a13, a23] = meshgrid(linspace(-range, range, N), linspace(-range, range, N));
a13sq = a13.^2;
a23sq = a23.^2;

%% Diverging blue-white-red colormap (stage-wide convention; see
% plot_topo_2node_dyad_subplot.m for the full rationale).
rdbu_anchors = [ ...
      5,  48,  97;   33, 102, 172;   67, 147, 195;  146, 197, 222; ...
    209, 229, 240;  247, 247, 247;  253, 219, 199;  244, 165, 130; ...
    214,  96,  77;  178,  24,  43;  103,   0,  31] / 255;
neutral_halfwidth = 0.05;
v = linspace(-1, 1, 256);
t = 0.5 + 0.5 * sign(v) .* max(abs(v) - neutral_halfwidth, 0) ...
                        ./ (1 - neutral_halfwidth);
cmap_div = interp1(linspace(0,1,size(rdbu_anchors,1)), rdbu_anchors, t);

rho_levels = -1:0.1:1;
rho_fill_levels = -1:0.02:1;

%% NO stability mask, NO grey invalid region, NO black boundary curve
% (topology-specific deviation, same as the combined confounder plotter): the
% confounder's A = [a,0,a13; 0,a,a23; 0,0,a] is upper triangular, eigenvalues
% all = a < 0, so the system is Hurwitz-stable everywhere in the plotted plane;
% and the denominators are strictly positive, so there is no boundary curve.

%% Output folder -- one subfolder per topology.
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fullfile(script_dir, '..', '..');
fig_dir    = fullfile(proj_root, 'figures', 'topo_3node_confounder_figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Fixed square canvas and fixed axes Position -- IDENTICAL for every panel
% (see plot_topo_2node_dyad_subplot.m for the full rationale on panel_size
% and ax_pos, and why exportgraphics is used).
panel_size = 360;               % figure width = height, in points (nominal)
ax_pos = [0.16 0.12 0.70 0.70]; % [left bottom width height], normalized

%% Per-panel: compute rho_* gradient, plot, export.
for k = 1:size(params,1)
    a        = params(k,1);
    sigma_W2 = params(k,2);
    sigma_W3 = params(k,3);

    % rho_21 for the confounder (verified against Appendix D
    % eq:d26_rho21_confounder / main text eq:rho21_confounder).
    numerator = a13 .* a23 * sigma_W3;
    denom1 = a13sq*sigma_W3 + 2*a*a*sigma_W1;
    denom2 = a23sq*sigma_W3 + 2*a*a*sigma_W2;
    rho_21 = numerator ./ sqrt(denom1 .* denom2);

    [drho_da13, drho_da23] = gradient(rho_21, a13(1,:), a23(:,1));

    % Standalone figure for this panel only -- 'Visible','off' avoids the
    % headless window-manager clamp.
    fig = figure('Color','w','Units','points', ...
        'Position',[100 100 panel_size panel_size],'Visible','off');
    ax = axes(fig,'Units','normalized','Position',ax_pos);

    contourf(ax, a13, a23, rho_21, rho_fill_levels, 'LineStyle', 'none')
    hold(ax,'on')
    contour(ax, a13, a23, rho_21, rho_levels(rho_levels ~= 0), ...
            'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.5)
    colormap(ax, cmap_div)
    clim(ax, [-1 1])

    % rho_21 = 0 contour, bright green dashed (stage-wide convention).
    contour(ax, a13, a23, rho_21, [0 0], '--', ...
            'LineColor', [0.00 0.62 0.05], 'LineWidth', 2.5)
    % (No black stability boundary: confounder is stable everywhere.)

    mag = sqrt(drho_da13.^2 + drho_da23.^2);
    U = drho_da13 ./ mag;
    V = drho_da23 ./ mag;
    step = 20;
    arrow_scale = 0.5;
    quiver(ax, a13(1:step:end,1:step:end), ...
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

    % Font sizes match the dyad reference: axis labels 22, title 18, ticks 16.
    xlabel(ax,'\textbf{$\mathbf{a_{13}}$}','Interpreter','latex','FontSize',22)
    ylabel(ax,'\textbf{$\mathbf{a_{23}}$}','Interpreter','latex','FontSize',22)

    % %+.2f (explicit sign) for a, matching the dyad convention, so titles have
    % uniform width. sigma values are always positive -> %.2f.
    title(ax, sprintf(['$a=%+.2f$, $\\sigma_{W_3}=%.2f$' newline ...
        '$\\sigma_{W_1}=%.2f$, $\\sigma_{W_2}=%.2f$'], ...
        a, sigma_W3, sigma_W1, sigma_W2), ...
        'Interpreter','latex','FontSize',18);

    hold(ax,'off')

    out_path = fullfile(fig_dir, [panel_labels{k} '.pdf']);
    exportgraphics(fig, out_path, 'ContentType', 'vector');
    fprintf('Wrote %s\n', out_path);
    close(fig);
end

%% Standalone colorbar (identical construction to the dyad reference; see
% plot_topo_2node_dyad_subplot.m for the no-dummy-image / AxisLocation /
% canvas-size rationale). Label is rho_21, same quantity the confounder plots.
cb_width  = 3 * panel_size * 1.1;       % 1188
cb_height = 160;                        % nominal only, autocropped down
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

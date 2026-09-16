%% Contour + vector-field plot of rho_21 for the 3-node mediator topology, as
% a function of a31, a23 (see figure caption for what each plotted element
% represents). Closed-form rho_21 derived in Appendix D (eq:rho21_mediator /
% main text eq:rho21_mediator_main).
%
% One standalone square PDF per panel, plus one standalone colorbar PDF,
% assembled in LaTeX per the authoritative guide
% (tasks/stage3_tasks_plot_networks/stage3_topo_x-node_figure_template_instructions.md).
% Structure copied from plot_topo_3node_confounder_subplot.m; only the
% topology-specific pieces differ (formula, axes, contour-line levels).
% Output: figures/topo_3node_mediator_figures/{B1,B2,B3,C1,C2,C3,colorbar}.pdf
clc; clear;

% Params for the mediator, panels B1-C3 (2x3 grid). columns: [a, sigma_W2,
% sigma_W3]; sigma_W1 = 1 throughout. Same six values as the collider figure
% (deliberately, to test the observation that the mediator's contour
% geometry looks like a "dampened version" of the confounder's -- though the
% confounder figure itself now uses different (a, sigma_W3) values, see that
% file). Row B fixes sigma_W2=sigma_W3=1 and varies a in {-0.3,-0.6,-1.0}
% across columns; row C varies sigma_W2 and/or sigma_W3 per column too.
params = [...
    -0.3, 1, 1;
    -0.6, 1, 1;
    -1.0, 1, 1;
    -0.3, 1, 3;
    -0.6, 3, 1;
    -1.0, 3, 3];
sigma_W1 = 1;
% Row A is reserved for the network diagram and matrices strip (drawn
% directly in LaTeX, not exported here), so the two MATLAB-generated data
% rows are B and C.
panel_labels = {'B1','B2','B3','C1','C2','C3'};

% param grid for a31 and a23 (the two mediation path coefficients: 1 -> 3 -> 2)
range = 1.5;
N = 150;
[a31, a23] = meshgrid(linspace(-range, range, N), linspace(-range, range, N));
a31sq = a31.^2;
a23sq = a23.^2;

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

rho_fill_levels = -1:0.02:1;
rho_levels = linspace(-1,1,31);  % 30 visible lines -- a denser 40-line
                                  % (0.05-step) version looked too busy

%% NO stability mask, NO grey invalid region, NO black boundary curve
% (topology-specific deviation, same as the combined mediator plotter):
% reordering the mediator's A = [a,0,0; 0,a,a23; a31,0,a] by node sequence
% (1,3,2) gives a lower-triangular matrix, so its eigenvalues ARE the diagonal
% entries (all = a < 0) -- Hurwitz-stable everywhere in the plotted plane. The
% denominator is a sum of non-negative terms, strictly positive whenever a!=0,
% so there is also no boundary curve.

%% Output folder -- one subfolder per topology.
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fullfile(script_dir, '..', '..');
fig_dir    = fullfile(proj_root, 'figures', 'topo_3node_mediator_figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Fixed square canvas and fixed axes Position -- IDENTICAL for every panel.
panel_size = 360;
ax_pos = [0.16 0.12 0.70 0.70];

%% Per-panel: compute rho_* gradient, plot, export.
for k = 1:size(params,1)
    a        = params(k,1);
    sigma_W2 = params(k,2);
    sigma_W3 = params(k,3);

    % rho_21 for the mediator (Appendix D eq:rho21_mediator / main text
    % eq:rho21_mediator_main).
    numerator = a31 .* a23 * sigma_W1;
    denom1 = 2*sigma_W1;
    denom2 = 4*a*a*a23sq*sigma_W3 + 8*a^4*sigma_W2 + 3*a23sq.*a31sq*sigma_W1;
    rho_21 = numerator ./ sqrt(denom1 .* denom2);

    [drho_da31, drho_da23] = gradient(rho_21, a31(1,:), a23(:,1));

    fig = figure('Color','w','Units','points', ...
        'Position',[100 100 panel_size panel_size],'Visible','off');
    ax = axes(fig,'Units','normalized','Position',ax_pos);

    contourf(ax, a31, a23, rho_21, rho_fill_levels, 'LineStyle', 'none')
    hold(ax,'on')
    contour(ax, a31, a23, rho_21, rho_levels(rho_levels ~= 0), ...
            'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.5)
    colormap(ax, cmap_div)
    clim(ax, [-1 1])

    % rho_21 = 0 contour, bright green dashed (stage-wide convention).
    contour(ax, a31, a23, rho_21, [0 0], '--', ...
            'LineColor', [0.00 0.62 0.05], 'LineWidth', 2.5)
    % (No black stability boundary: mediator is stable everywhere.)

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

    % Font sizes match the dyad/confounder/collider reference: axis labels 22,
    % title 18, ticks 16.
    xlabel(ax,'\textbf{$\mathbf{a_{31}}$}','Interpreter','latex','FontSize',22)
    ylabel(ax,'\textbf{$\mathbf{a_{23}}$}','Interpreter','latex','FontSize',22)

    % %+.2f (explicit sign) for a, matching the confounder/collider convention.
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

%% Standalone colorbar (identical construction to the dyad/confounder/collider
% reference). Label is rho_21, the quantity the mediator plots.
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

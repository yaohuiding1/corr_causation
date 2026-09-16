%% Contour + vector-field plot of rho_21 for the 3-node two-way mediator
% topology, as a function of a31, a23 (see figure caption for plotted
% elements). Closed-form rho_21 derived in Appendix D (eq:rho21_mediator2way
% / main text eq:rho21_med2way_main).
%
% One standalone square PDF per panel, plus one standalone colorbar PDF,
% assembled in LaTeX per the authoritative guide
% (tasks/stage3_tasks_plot_networks/stage3_topo_x-node_figure_template_instructions.md).
% Structure copied from plot_topo_2node_dyad_subplot.m, since this topology
% (unlike confounder/collider/mediator) is cyclic and genuinely needs a
% stability mask.
% Output: figures/topo_3node_mediator_2way_figures/{B1,B2,B3,C1,C2,C3,colorbar}.pdf
%
% The six panels don't split into a clean 3+3 row contrast (panels 1-4 share
% a=-0.5; panels 5-6 use a=-1.0, sigma_W3=3); preserved as-is from the
% original combined figure, relabeled B1-C3 in sequence.
%
% No sign(K) factor is applied: the main text (corr_causation_v4.tex,
% Two-Way Mediator subsection) proves |K| = 16a^2(s-a^2)(s-4a^2), with
% s=a13*a31+a23*a32, is positive throughout the stability region, so
% sign(K)=+1 identically there.
clc; clear;

% Six panels: a13, a32, a, sigma_W1, sigma_W2, sigma_W3 (unchanged from the
% combined script). Panels A1-A4-equivalent (1-4) share a=-0.5, sigma_W=(1,1,1)
% and rotate (a13,a32) through sign quadrants (starting at (0,0), which
% recovers the one-way mediator exactly); panels 5-6 (B2,B3) repeat panels
% 3-4's (a13,a32) pairs with a=-1.0, sigma_W3=3 instead, isolating how node 3's
% own dynamics reshape the geometry for a fixed pair of return connections.
params = [...
%    a13   a32     a   sW1  sW2  sW3
     0.0,  0.0,  -0.5,   1,   1,   1;
     1.0,  1.0,  -0.5,   1,   1,   1;
    -1.0,  1.0,  -0.5,   1,   1,   1;
    -1.0, -1.0,  -0.5,   1,   1,   1;
    -1.0,  1.0,  -1.0,   1,   1,   3;
    -1.0, -1.0,  -1.0,   1,   1,   3];
% Row A is reserved for the network diagram and matrices strip (drawn
% directly in LaTeX, not exported here), so the two MATLAB-generated data
% rows are B and C.
panel_labels = {'B1','B2','B3','C1','C2','C3'};

% param grid for a31 and a23 (the swept axes; a13, a32 are fixed per panel).
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
rho_levels = -1:0.1:1;
rho_fill_levels = -1:0.02:1;

%% Output folder -- one subfolder per topology.
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fullfile(script_dir, '..', '..');
fig_dir    = fullfile(proj_root, 'figures', 'topo_3node_mediator_2way_figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Fixed square canvas and fixed axes Position -- IDENTICAL for every panel.
panel_size = 360;
ax_pos = [0.16 0.12 0.70 0.70];

%% Per-panel: compute rho_* gradient, plot, export.
for k = 1:size(params,1)
    a13 = params(k,1);
    a32 = params(k,2);
    a   = params(k,3);
    sigma_W1 = params(k,4);
    sigma_W2 = params(k,5);
    sigma_W3 = params(k,6);

    % rho_21 for the two-way mediator (verified against Appendix D
    % eq:rho21_mediator2way / main text eq:rho21_med2way_main). sign(K)
    % dropped -- see file header note.
    N21 = 2*a^2*sigma_W1*(a23.*a31) + a13*sigma_W1*(a23.*a31sq) ...
        - 2*a32*sigma_W1*(a23sq.*a31) ...
        + 2*a^2*a13*a32*sigma_W2 - 2*a13^2*a32*sigma_W2*a31 ...
        + a13*a32^2*sigma_W2*a23 ...
        + 4*a^2*a13*sigma_W3*a23 - a13^2*sigma_W3*(a23.*a31) ...
        - a13*a32*sigma_W3*a23sq;

    D11 = 8*a^4*sigma_W1 - 6*a^2*a13*sigma_W1*a31 - 10*a^2*a32*sigma_W1*a23 ...
        + a13^2*sigma_W1*a31sq + 2*a32^2*sigma_W1*a23sq ...
        + 3*a13^2*a32^2*sigma_W2 ...
        + 4*a^2*a13^2*sigma_W3 - a13^3*sigma_W3*a31 - a13^2*a32*sigma_W3*a23;

    D22 = 3*sigma_W1*(a23sq.*a31sq) ...
        + 8*a^4*sigma_W2 - 10*a^2*a13*sigma_W2*a31 - 6*a^2*a32*sigma_W2*a23 ...
        + 2*a13^2*sigma_W2*a31sq + a32^2*sigma_W2*a23sq ...
        + 4*a^2*sigma_W3*a23sq - a32*sigma_W3*a23cube - a13*sigma_W3*(a23sq.*a31);

    rho_21 = N21 ./ sqrt(D11 .* D22);

    % Stability mask: a13*a31 + a23*a32 < a^2 (Hurwitz condition, a<0 always
    % true here).
    valid = (a13*a31 + a23*a32) < a^2;
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
    % Stability boundary, black solid.
    contour(ax, a31, a23, a13*a31 + a23*a32, [a^2 a^2], 'k-', 'LineWidth', 2.5)

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

    % Font sizes match the dyad/confounder/collider/mediator reference: axis
    % labels 22, title 18, ticks 16.
    xlabel(ax,'\textbf{$\mathbf{a_{31}}$}','Interpreter','latex','FontSize',22)
    ylabel(ax,'\textbf{$\mathbf{a_{23}}$}','Interpreter','latex','FontSize',22)

    % Title shows ALL four varying quantities (a13, a32, a, sigma_W3) -- the
    % combined script's title showed only a13/a32, relying on the caption to
    % explain that a/sigma_W3 also change for panels 5-6; made explicit here
    % to match this template's per-panel-title convention (every varying
    % parameter shown in the title itself, not just the caption).
    title(ax, sprintf(['$a_{13}=%+.2f$, $a_{32}=%+.2f$' newline ...
        '$a=%+.2f$, $\\sigma_{W_3}=%.2f$'], ...
        a13, a32, a, sigma_W3), ...
        'Interpreter','latex','FontSize',18);

    hold(ax,'off')

    out_path = fullfile(fig_dir, [panel_labels{k} '.pdf']);
    exportgraphics(fig, out_path, 'ContentType', 'vector');
    fprintf('Wrote %s\n', out_path);
    close(fig);
end

%% Standalone colorbar (identical construction to the dyad/confounder/
% collider/mediator reference). Label is rho_21, the quantity plotted here.
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

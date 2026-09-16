%% Contour + vector-field plot of rho_21 for the 2-node dyad topology, as a
% function of the causal parameters a12, a21 (see figure caption for what
% each plotted element represents). Closed-form rho_21 derived in Appendix C
% (appendices/appendix_c_dyadic_interaction.tex).
%
% Exports each panel as its own standalone square PDF, plus one standalone
% colorbar PDF, so LaTeX assembles the Figure 2 grid from independent,
% uniformly-sized pieces (rather than one combined tiledlayout image).
% Output: figures/topo_2node_dyad_figures/{B1,B2,B3,C1,C2,C3,colorbar}.pdf
%
% LAYOUT: each panel is its own figure with NO tiledlayout, and uses an
% EXPLICIT, IDENTICAL axes Position (normalized units) inside an identical
% square figure canvas. A plain axes with a fixed Position renders at a fixed
% size regardless of content, so all six exports are guaranteed geometrically
% identical -- not just visually close.
clc; clear;

% Six (a11, a22, sigma_w1, sigma_w2) combinations, panels B1-C3. Row 3 /
% column 3 intentionally swaps the sign of a11/a22 relative to the other
% rows in that column.
params = [...
    -0.5, -0.5, 1, 1;
    -1.5, -0.5, 1, 1;
     0.3, -0.8, 1, 1.0;
    -0.5, -0.5, 3, 1;
    -1.5, -0.5, 3, 1;
    -0.8,  0.3, 3, 1];
% Row A is reserved for the network diagram and matrices strip (drawn
% directly in LaTeX, not exported here), so the two MATLAB-generated data
% rows are B and C.
panel_labels = {'B1','B2','B3','C1','C2','C3'};

% param grid for a12 and a21
range = 1.5;
N = 150;
[a12, a21] = meshgrid(linspace(-range, range, N), linspace(-range, range, N));
a12sq = a12.^2;
a21sq = a21.^2;
a12a21 = a12.*a21;

%% Diverging blue-white-red colormap, with a widened flat/neutral band near
% zero (rho ~= 0 maps to white over a wider range than a plain linear map).
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

%% Output folder -- one subfolder per topology, individual panels + colorbar
% inside.
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fullfile(script_dir, '..', '..');
fig_dir    = fullfile(proj_root, 'figures', 'topo_2node_dyad_figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Fixed square canvas and fixed axes Position -- IDENTICAL for every panel.
% Figure: 360x360 points (nominal, before autocrop). Axes box: 70% of the
% canvas (square), leaving room for a 2-line title above and labels/ticks
% below and left. Units are 'points' for an exact mapping to PDF points.
%
% Uses exportgraphics (sharper text than print('-painters'), matching the
% combined figure's pipeline) over exact identical sizes -- it autocrops to
% ink, so panel PDFs vary by a few points with title-text width; acceptable
% since panels don't need to be pixel-identical. panel_size=360 keeps the
% autocropped ink width close to the combined figure's per-panel width
% (~340pt), so the same absolute FontSize looks similarly sized in both.
panel_size = 360;               % figure width = height, in points (nominal,
                                 % before exportgraphics's autocrop)
ax_pos = [0.16 0.12 0.70 0.70]; % [left bottom width height], normalized;
                                 % width==height -> axes box is square since
                                 % the figure itself is square.

%% Per-panel: compute rho_* gradient, plot, export.
for k = 1:size(params,1)
    a11 = params(k,1);
    a22 = params(k,2);
    sigma_w1 = params(k,3);
    sigma_w2 = params(k,4);

    % Closed-form rho_21 for the 2-node dyad.
    numerator = -(a22.*a21*sigma_w1 + a11.*a12*sigma_w2);
    bracket_term1 = a22*(a11+a22) - a12a21;
    bracket_term2 = a11*(a11+a22) - a12a21;
    denom1 = bracket_term1*sigma_w1 + a12sq*sigma_w2;
    denom2 = bracket_term2*sigma_w2 + a21sq*sigma_w1;
    rho_21 = numerator ./ sqrt((denom1 .* denom2));

    % Hurwitz stability mask (2x2 case: trace<0 and det>0).
    trace_ok = (a11 + a22) < 0;
    det_ok   = a12a21 < a11*a22;
    valid    = trace_ok & det_ok;

    rho_21_masked = rho_21;
    rho_21_masked(~valid) = NaN;

    [drho_da12, drho_da21] = gradient(rho_21_masked, a12(1,:), a21(:,1));
    drho_da12(~valid) = NaN;
    drho_da21(~valid) = NaN;

    % Standalone figure for this panel only -- 'Visible','off' avoids the
    % headless window-manager clamp.
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

    contour(ax, a12, a21, rho_21_masked, [0 0], '--', ...
            'LineColor', [0.00 0.62 0.05], 'LineWidth', 2.5)
    contour(ax, a12, a21, a12a21, [a11*a22 a11*a22], 'k-', 'LineWidth', 2.5)

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

    % FontSize must be set before xlabel/ylabel/title, or they inherit the
    % axes' default size instead of this one.
    set(ax,'FontSize',16)

    axis_ticks = -range:0.5:range;
    set(ax, 'XTick', axis_ticks, 'YTick', axis_ticks)
    ax.XAxis.TickLabelRotation = 0;
    ax.YAxis.TickLabelRotation = 0;

    % Axis-label FontSize (22) matches the visual weight of the 18pt title
    % once LaTeX-rendered; tick FontSize stays at 16. Sizes run larger than
    % usual since panels are shrunk substantially when embedded in the grid.
    xlabel(ax,'\textbf{$\mathbf{a_{12}}$}','Interpreter','latex','FontSize',22)
    ylabel(ax,'\textbf{$\mathbf{a_{21}}$}','Interpreter','latex','FontSize',22)

    % %+.2f for a11/a22 (can be negative) keeps title width -- and thus
    % autocrop -- consistent across panels; sigma_w1/sigma_w2 (always
    % positive) don't need it. sigma_W_ is capitalized to match Sigma_W,
    % same convention as sigma_11 for Sigma_Y.
    title(ax, sprintf(['$a_{11}=%+.2f$, $a_{22}=%+.2f$' newline ...
        '$\\sigma_{W_1}=%.2f$, $\\sigma_{W_2}=%.2f$'], ...
        a11, a22, sigma_w1, sigma_w2), ...
        'Interpreter','latex','FontSize',18);

    hold(ax,'off')

    out_path = fullfile(fig_dir, [panel_labels{k} '.pdf']);
    exportgraphics(fig, out_path, 'ContentType', 'vector');
    fprintf('Wrote %s\n', out_path);
    close(fig);
end

%% Standalone colorbar.
%
% colorbar() only needs its parent axes' Colormap/CLim, not a plotted image,
% so the axes stays empty and invisible -- a dummy image here would stretch
% to fill the canvas (no fixed aspect ratio) and produce an unwanted
% tall/stretched rectangle behind the bar.
%
% CANVAS SIZE: width = 3 * panel_size, plus 10% margin, to stretch the bar
% wider relative to the fixed-point-size tick/label text once autocropped
% and shrunk to the embedded row width. Height is just a nominal upper
% bound -- with no dummy image, exportgraphics autocrops to just the bar +
% ticks + label.
cb_width  = 3 * panel_size * 1.1;       % 1188
cb_height = 160;                        % nominal only, autocropped down
fig = figure('Color','w','Units','points','Position',[100 100 cb_width cb_height], ...
    'Visible','off');
ax = axes(fig,'Visible','off');
colormap(ax, cmap_div)
clim(ax, [-1 1])
cb = colorbar(ax,'Location','south');
% AxisLocation='out' (default is 'in'): puts ticks/label on the outward-
% facing side of the bar. With the default 'in', ticks would face into the
% axes and appear to sit on top of the bar.
cb.AxisLocation = 'out';
% Bar thickness (Position height) kept small (0.20) since exportgraphics
% autocrops the rest to ink anyway -- this is the actual height lever.
cb.Position = [0.08 0.35 0.84 0.20];
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathbf{\rho_{21}}$'; % already bold via \mathbf --
                                          % Interpreter='latex' ignores
                                          % FontWeight on the Text object
cb.Label.FontWeight = 'bold';
cb.FontSize = 18;
cb.Label.FontSize = 22;

cb_path = fullfile(fig_dir, 'colorbar.pdf');
exportgraphics(fig, cb_path, 'ContentType', 'vector');
fprintf('Wrote %s\n', cb_path);
close(fig);

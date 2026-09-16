%% Contour + vector-field plot of sigma_33 (raw variance of node 3) for the
% 3-node coupled-drivers topology, as a function of the reciprocal-dyad
% parameters a12, a21 (see figure caption for plotted elements). Plots
% sigma_33 (Appendix D, eq:sigma33_coupled), not rho_21: rho_21 here equals
% the dyadic formula exactly (eq:rho21_coupled_main = eq:rho21_dyad at
% a11=a22=a), so it would be redundant with Figure 2.
%
% One standalone square PDF per panel, plus a colorbar PDF, assembled in
% LaTeX per the authoritative guide
% (tasks/stage3_tasks_plot_networks/stage3_topo_x-node_figure_template_instructions.md).
% Output: figures/topo_3node_coupled_drivers_figures/{B1,B2,B3,C1,C2,C3,colorbar}.pdf
%
% Sweeps (a12,a21); (a31,a32,sigma_W3) vary per panel. Verified
% symbolically: (1) a stability mask IS needed, since the boundary
% hyperbola a12*a21=a^2 cuts through the swept plane (~68% stable at
% a=-0.5). (2) sigma_33 depends on a31,a32 only via a31^2, a32^2, a31*a32,
% so it's invariant under sign flips of both; the six (a31,a32) sets below
% include a rotation pair, (0.3,0.3) vs (0.3,-0.3), that displays this.
%
% Uses a sequential viridis colormap (not diverging, no neutral plateau, no
% zero-contour) since sigma_33 is a non-negative, unbounded variance, not a
% bounded correlation. Everything else carries over from the rho-plotting
% template (shared clim, quiver field, ticks, fonts, panel geometry, export
% mechanics).
clc; clear;

a = -0.5;
sigma_W1 = 1; sigma_W2 = 1;

% Six (a31,a32) sets, sigma_W3 fixed at 1: the coefficient of sigma_W3 in
% sigma_33 simplifies exactly to -1/(2a) -- a constant, independent of
% a12,a21,a31,a32 (verified symbolically) -- so sigma_W3 only translates the
% surface, never reshapes it. A sigma_W3-contrast row (as used for the other
% topologies) would be informationally empty here, so all six panels instead
% sweep distinct (a31,a32) values.
%
% a31,a32 do NOT affect stability (they are not swept, and never entered |K|),
% so any values are admissible. The six below span the three structural
% degrees of freedom -- sigma_33 depends on a31,a32 only through a31^2, a32^2
% and the cross term a31*a32:
params = [...
%    a31    a32   sigma_W3
     0.3,   0.3,   1;   % B1: equal magnitude, same sign      (cross +0.09)
     0.3,  -0.3,   1;   % B2: equal magnitude, opposite sign  (cross -0.09)
                        %     = B1 rotated 180 deg -- shows the sign structure
     0.3,   0.5,   1;   % B3: unequal magnitude, same sign    (cross +0.15)
     0.0,   0.3,   1;   % C1: single driver only, cross term vanishes
     0.5,   0.5,   1;   % C2: equal magnitude, same sign, stronger (cross +0.25)
     0.6,  -0.2,   1];  % C3: strongly unequal, opposite sign (cross -0.12)
% Row A is reserved for the network diagram and matrices strip (drawn
% directly in LaTeX, not exported here), so the two MATLAB-generated data
% rows are B and C.
panel_labels = {'B1','B2','B3','C1','C2','C3'};

% param grid for a12 and a21 (the reciprocal dyad between the two drivers).
% Range [-1.5,1.5], matching the stage-wide convention.
range = 1.5;
N = 150;
[a12, a21] = meshgrid(linspace(-range, range, N), linspace(-range, range, N));
a12a21 = a12 .* a21;
valid  = a12a21 < a^2;          % Hurwitz condition (a<0 fixed above)
invalid_grey = [0.80 0.80 0.80];

%% Sequential viridis-like colormap (not built into this MATLAB install;
% constructed from the standard published viridis anchor RGB values, same
% anchors+interp1 pattern as the diverging colormap used elsewhere).
viridis_anchors = [ ...
     68,   1,  84;   72,  40, 120;   62,  74, 137;   49, 104, 142; ...
     38, 130, 142;   31, 158, 137;   53, 183, 121;  109, 205,  89; ...
    253, 231,  37] / 255;
cmap_seq = interp1(linspace(0,1,size(viridis_anchors,1)), viridis_anchors, linspace(0,1,256));

%% sigma_33, as a function of the swept (a12,a21) for one (a31,a32,sigma_W3).
% Appendix D eq:sigma33_coupled -- verified term-for-term from the appendix.
sigma33_fn = @(a31,a32,sw3) ...
    -(  (16*a^3*a31^2 - 24*a^2*a21*a31*a32 + 12*a*a21.^2*a32^2 - 4*a*a12a21*a31^2) * sigma_W1 ...
      + (16*a^3*a32^2 - 24*a^2*a12*a31*a32 + 12*a*a12.^2*a31^2 - 4*a*a12a21*a32^2) * sigma_W2 ...
      + (32*a^5 - 40*a^3*a12a21 + 8*a*a12a21.^2) * sw3 ) ...
    ./ (16*a^2*(a^2-a12a21).*(4*a^2-a12a21));

%% Shared LOG color range across all six panels.
% sigma_33 DIVERGES at the stability boundary a12*a21 -> a^2 (it is a
% steady-state variance, and the steady state ceases to exist there), so a
% max-based clim is meaningless: the observed max is purely an artefact of how
% close a grid point lands to the boundary hyperbola, and grows without bound
% with N (sweep: N=150 -> 9.6e3, N=300 -> 1.13e5, N=600 -> 1.60e5).
%
% Upper limit: a fixed cap of sigma_33 = 3.5. Without a cap, contour lines
% bunch tightly near the stability boundary (sigma_33 diverges there), so
% most of the 13 log-spaced levels get spent resolving a thin high-value
% sliver, leaving the broad low/mid-value plateau (most of each panel's
% area) with few or no lines; lowering the cap removes those high,
% tightly-packed levels (diagnostic sweep, panel B1: cap 5 -> 74.9% of area
% in the bottom 3 of 13 bands, cap 3.5 -> 67.7%, cap 3 -> 63.2%). Raising the
% lower limit instead was rejected: the low values sit on a wide flat
% plateau, so nudging the floor from 1.05 to 1.2 would flatten ~22% of each
% panel to one dead colour -- a worse trade than the top-end saturation
% increase from 5->3.5. Lower limit is just the pooled minimum (the floor is
% set by sigma_W3).
pool = [];
for k = 1:size(params,1)
    S = sigma33_fn(params(k,1), params(k,2), params(k,3));
    S(~valid) = NaN;
    pool = [pool; S(:)]; %#ok<AGROW>
end
pool = pool(isfinite(pool) & pool > 0);
sigma33_cap = 3.5;
clim_lo = log10(min(pool));
clim_hi = log10(sigma33_cap);
fprintf('sigma_33 pooled over stable region: min %.3f  median %.3f  p90 %.3f\n', ...
        min(pool), median(pool), prctile(pool,90));
fprintf('cap = %g (log10 %.4f); %.2f%% of stable plane saturates\n', ...
        sigma33_cap, clim_hi, 100*mean(pool > sigma33_cap));
fill_levels = linspace(clim_lo, clim_hi, 100);
line_levels = linspace(clim_lo, clim_hi, 13);

%% Output folder -- one subfolder per topology.
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fullfile(script_dir, '..', '..');
fig_dir    = fullfile(proj_root, 'figures', 'topo_3node_coupled_drivers_figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Fixed square canvas and fixed axes Position -- IDENTICAL for every panel.
panel_size = 360;
ax_pos = [0.16 0.12 0.70 0.70];

%% Per-panel: compute sigma_33 gradient, plot, export.
for k = 1:size(params,1)
    a31      = params(k,1);
    a32      = params(k,2);
    sigma_W3 = params(k,3);

    sigma33 = sigma33_fn(a31, a32, sigma_W3);
    sigma33(~valid) = NaN;
    log_sigma33 = log10(sigma33);

    % Gradient direction is identical whether computed on sigma_33 or
    % log10(sigma_33) (log is monotonic increasing on a strictly positive
    % quantity, so grad(log f) = grad(f)/f, a positive scalar multiple --
    % same direction). Computed on the raw quantity for simplicity.
    [dsig_da12, dsig_da21] = gradient(sigma33, a12(1,:), a21(:,1));
    dsig_da12(~valid) = NaN;
    dsig_da21(~valid) = NaN;

    fig = figure('Color','w','Units','points', ...
        'Position',[100 100 panel_size panel_size],'Visible','off');
    ax = axes(fig,'Units','normalized','Position',ax_pos);

    contourf(ax, a12, a21, log_sigma33, fill_levels, 'LineStyle', 'none')
    hold(ax,'on')
    contour(ax, a12, a21, log_sigma33, line_levels, ...
            'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.5)
    colormap(ax, cmap_seq)
    set(ax, 'Color', invalid_grey)
    clim(ax, [clim_lo clim_hi])

    % Stability boundary, black solid (stage-wide convention).
    contour(ax, a12, a21, a12a21, [a*a a*a], 'k-', 'LineWidth', 2.5)

    mag = sqrt(dsig_da12.^2 + dsig_da21.^2);
    U = dsig_da12 ./ mag;
    V = dsig_da21 ./ mag;
    step = 20;
    arrow_scale = 0.5;
    quiver(ax, a12(1:step:end,1:step:end), ...
        a21(1:step:end,1:step:end), ...
        U(1:step:end,1:step:end), ...
        V(1:step:end,1:step:end), ...
        arrow_scale,'w', 'LineWidth', 1)  % white arrows: viridis is dark at
                                           % low values, black would vanish

    axis(ax,'square')

    % set FontSize before xlabel/ylabel/title (ordering pitfall).
    set(ax,'FontSize',16)

    axis_ticks = -range:0.5:range;
    set(ax, 'XTick', axis_ticks, 'YTick', axis_ticks)
    ax.XAxis.TickLabelRotation = 0;
    ax.YAxis.TickLabelRotation = 0;

    % Font sizes match the rho-figure reference: axis labels 22, title 18,
    % ticks 16.
    xlabel(ax,'\textbf{$\mathbf{a_{12}}$}','Interpreter','latex','FontSize',22)
    ylabel(ax,'\textbf{$\mathbf{a_{21}}$}','Interpreter','latex','FontSize',22)

    % %+.2f (explicit sign) for a31/a32, matching the stage-wide sign
    % convention. sigma_W3 is fixed across all six panels, but the second
    % title line is kept anyway so panel geometry matches the other
    % topologies' two-line titles (row heights are matched in the assembled
    % LaTeX figure).
    title(ax, sprintf(['$a_{31}=%+.2f$, $a_{32}=%+.2f$' newline ...
        '$\\sigma_{W_3}=%.2f$'], a31, a32, sigma_W3), ...
        'Interpreter','latex','FontSize',18);

    hold(ax,'off')

    out_path = fullfile(fig_dir, [panel_labels{k} '.pdf']);
    exportgraphics(fig, out_path, 'ContentType', 'vector');
    fprintf('Wrote %s\n', out_path);
    close(fig);
end

%% Standalone colorbar (same mechanics as the rho-figure reference: no dummy
% image, AxisLocation='out', position/font sizes). Label is sigma_33 (not
% rho_21); colormap is the sequential viridis-like map, LOG scale. Ticks are
% placed at real sigma_33 values (not raw log10 numbers) via explicit
% Ticks/TickLabels, chosen as round numbers inside the plotted range.
cb_width  = 3 * panel_size * 1.1;
cb_height = 160;
fig = figure('Color','w','Units','points','Position',[100 100 cb_width cb_height], ...
    'Visible','off');
ax = axes(fig,'Visible','off');
colormap(ax, cmap_seq)
clim(ax, [clim_lo clim_hi])
cb = colorbar(ax,'Location','south');
cb.AxisLocation = 'out';
cb.Position = [0.08 0.35 0.84 0.20];
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathbf{\sigma_{33}}$';
cb.Label.FontWeight = 'bold';
cb.FontSize = 18;
cb.Label.FontSize = 22;
cand = [0.5 1 2 3 5 10 20 30 50 100 200 300 500 1000];
tick_vals = cand(log10(cand) >= clim_lo & log10(cand) <= clim_hi);
cb.Ticks = log10(tick_vals);
cb.TickLabels = arrayfun(@(v) sprintf('%g', v), tick_vals, 'UniformOutput', false);

cb_path = fullfile(fig_dir, 'colorbar.pdf');
exportgraphics(fig, cb_path, 'ContentType', 'vector');
fprintf('Wrote %s\n', cb_path);
close(fig);

# Appendix D Figure Layout Spacing Parameters

All figures use a **1×4 horizontal layout**:
`[Diagram] gap0 [A] gap1 [Σ_Y] gap2 [Σ_W]`

All widths are fractions of `\textwidth`.

---

## TikZ Diagram Parameters (shared across all files)

| Parameter       | Value       | Description                        |
|-----------------|-------------|------------------------------------|
| `\a`            | `1.2 cm`    | Triangle size (node spacing)       |
| `\nodeSize`     | `0.7`       | Node circle diameter (cm)          |
| `\r`            | `0.5*\nodeSize cm` = `0.35 cm` | Self-loop radius  |
| Matrix font     | `{\large}`  | Font size for all three matrices   |

---

<!--
## 2-Node Topology

| File                  | Topology        | Diagram | gap0   | **A** | gap1   | **Σ_Y** | gap2   | **Σ_W** | Row total |
|-----------------------|-----------------|---------|--------|-------|--------|---------|--------|---------|-----------|
| `topo_2node_dyad`   | Reciprocal dyad | 0.25    | 0.02   | 0.18  | 0.02   | 0.18    | 0.04   | 0.18    | **0.87**  |
-->

---

## 3-Node Topologies

All 3-node files share identical spacing parameters:

| Parameter | Value          |
|-----------|----------------|
| Diagram   | `0.23\textwidth` |
| gap0 (diagram → A)   | `0\textwidth`    |
| A matrix  | `0.20\textwidth` |
| gap1 (A → Σ_Y)       | `0.06\textwidth` |
| Σ_Y matrix | `0.20\textwidth` |
| gap2 (Σ_Y → Σ_W)    | `0.08\textwidth` |
| Σ_W matrix | `0.20\textwidth` |
| **Row total** | **0.97\textwidth** |

Files using these parameters:

| File                        | Topology                      |
|-----------------------------|-------------------------------|
| `topo_3node_mediator`       | Indirect mediator (1→3→2)     |
| `topo_3node_mediator_2way`  | Bidirectional mediator        |
| `topo_3node_confounder`     | Confounder (common cause)     |
| `topo_3node_moderator`      | Reciprocal + common driver    |
| `topo_3node_collider`       | Collider (1→3←2)              |
| `topo_3node_coupled_drivers`| Coupled drivers               |
| `topo_3node_cycle`          | Directed cycle                |
| `topo_3node_full`           | Fully connected 3×3           |

---

## Notes

- **Interword space suppression**: in all 3-node files, `\end{minipage}` and
  `\hspace{...}` lines end with `%` to prevent TeX from inserting interword
  spaces between minipages. Without this, the row total would effectively
  exceed `\textwidth` and cause the elements to wrap.
- **2-node file** does not need `%` suppression (row total = 0.87\textwidth,
  well within budget even with interword spaces).
- **`topo_3node_full`** contains a second (legacy) figure layout using
  `{0.4\textwidth}` minipages; the spacing table above refers to the
  1×4 layout (Layout B) only.
- **Self-loop macro** `\selfloop{node}{pos}` is defined in
  `fig1_network_diagram_macros.tex`. Positions: 1 = left, 2 = right, 3 = top.
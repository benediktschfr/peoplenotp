# Architectural Specification: `peoplenotp` Clinical Certainty Suite

## 1. Overview & Design Philosophy
Build a cohesive, modular suite of functions and S3 methods for the `peoplenotp` R package to evaluate and visualize Bayesian clinical significance at the **population (group)** level, the **individual (person)** level, and across **sensitivity grids**.

**Core Workflow Principle:**
- **Group Level (`group_certainty()`):** Operates directly on a fitted `brmsfit` object (or posterior draws) to summarize fixed effects and inter-individual random-effects heterogeneity.
- **Individual Level (`individual_certainty()`):** Operates strictly on a user-prepared tidy data frame of posterior draws. Generates an S3 object containing an inverted Makowski Evidence Matrix.
- **Sensitivity Analysis:** Utilizes the existing `get_clinical_sensitivity_grid()` to produce a tibble of cohort success rates across a grid of ROPEs and certainty thresholds, visualized via `plot_clinical_sensitivity()`.

---

## 2. Function 1: `group_certainty()`

### Signature
```r
group_certainty(
  model,
  rope,
  better_is = c("lower", "higher"),
  parameters = NULL,
  ci = 0.95,
  ...
)
```

### Arguments
- `model`: A `brmsfit` object.
- `rope`: Numeric vector of length 2 specifying the Region of Practical Equivalence, e.g., `c(-2.5, 2.5)`.
- `better_is`: Character string indicating the clinical polarity of the outcome. 
  - `"lower"` (default): Negative change/lower values indicate improvement (e.g., PHQ-9, Burnout).
  - `"higher"`: Positive change/higher values indicate improvement (e.g., active minutes, Morning Energy).
- `parameters`: Optional character vector or regex pattern to filter target parameters.
- `ci`: Credible interval level for parameter point-and-interval summaries (default `0.95`). 

### Functional Requirements
1. **Directional Full Posterior Partitioning:**
   - Partition **100% of the full posterior distribution** for each fixed effect into three mutually exclusive directional columns (summing to 100%) based on the `better_is` argument:
     - If `better_is == "lower"`: 
       - `improvement_pct` = $P(\theta < \text{ROPE}_{\min})$
       - `rope_pct` = $P(\text{ROPE}_{\min} \le \theta \le \text{ROPE}_{\max})$
       - `deterioration_pct` = $P(\theta > \text{ROPE}_{\max})$
     - If `better_is == "higher"`:
       - `improvement_pct` = $P(\theta > \text{ROPE}_{\max})$
       - `rope_pct` = $P(\text{ROPE}_{\min} \le \theta \le \text{ROPE}_{\max})$
       - `deterioration_pct` = $P(\theta < \text{ROPE}_{\min})$
   - Also compute the Probability of Direction (`pd`), Posterior Median, and `95% CrI` (`CI_low`, `CI_high`).
2. **Random Effects (Individual Heterogeneity):**
   - Extract the group-level standard deviations (`sd__*`) and correlations (`cor__*`) from the `brmsfit` object (Posterior Median and 95% CrI).
3. **Return Object:** Return an S3 object of class `c("pnp_group_certainty", "data.frame")` and implement `print.pnp_group_certainty()` using `insight::export_table()`.

---

## 3. Function 2: `individual_certainty()`

### Signature
```r
individual_certainty(
  data,
  id,
  value,
  mid,
  better_is = c("lower", "higher"),
  draw_col = .draw,
  ci = 0.95
)
```

### Arguments
- `data`: A tidy data frame/tibble containing posterior draws per individual.
- `id`: Unquoted column name uniquely identifying individuals.
- `value`: Unquoted column name containing the individual target estimate per draw (e.g., `.epred`, `net_change`).
- `mid`: Numeric scalar (e.g., `2.5`) or numeric vector of positive values (e.g., `seq(1.5, 3.5, by = 0.5)`) defining the Minimal Important Difference threshold(s).
- `better_is`: Character string indicating clinical polarity (`"lower"` or `"higher"`).
- `draw_col`: Unquoted column name identifying the MCMC draw index (default `.draw`).
- `ci`: Credible interval width for individual summaries.

### Internal Computations (Executed across every value in `mid`)
1. **Individual-Level Posterior Summaries (`$individual_summary`):**
   - Calculate `median`, `ci_low`, `ci_high`, and `pd`.
   - Calculate `rope_pct`, `improvement_pct`, and `deterioration_pct` exactly as defined in `group_certainty()` based on `better_is`.
   - `directed_certainty`: $\max(\text{improvement\_pct}, \text{deterioration\_pct})$. 
   - `makowski_tier`: Ordered factor based on ROPE probability:
     - `"Negligible (>99% ROPE)"`
     - `"Prob. Negligible (>97.5%)"`
     - `"Undecided"`
     - `"Prob. Significant (<2.5%)"`
     - `"Significant (<1% ROPE)"`

2. **Inverted Evidence Matrix (`$matrix_data`):**
   - Tabulate counts ($n$) and proportions ($\%$) crossing `makowski_tier` (Rows) by `mid` thresholds (Columns).

### Return Structure
Return an S3 list object of class `"pnp_individual_certainty"` containing `$individual_summary`, `$matrix_data`, and `$meta`.

---

## 4. Print & Summary Methods (`print.pnp_individual_certainty`)

### Output Formatting
The print method exclusively outputs the inverted Clinical Significance Evidence Matrix using `insight::export_table()`. This prevents console clutter and fits perfectly into book margins. 

By default (`format = "both"`), the matrix combines counts and proportions into a single `n (%)` string per cell. 

**Expected Output (`format = "both"`):**
```text
# Clinical Significance Evidence Matrix

Evidence Tier               | MID = 2.0  | MID = 2.5  | MID = 3.0
------------------------------------------------------------------
Negligible (>99% ROPE)      |  2 (6.7%)  |  4 (13.3%) |  7 (23.3%)
Prob. Negligible (>97.5%)   |  4 (13.3%) |  5 (16.7%) |  6 (20.0%)
Undecided                   |  6 (20.0%) |  7 (23.3%) |  6 (20.0%)
Prob. Significant (<2.5%)   |  5 (16.7%) |  4 (13.3%) |  4 (13.3%)
Significant (<1% ROPE)      | 13 (43.3%) | 10 (33.3%) |  7 (23.3%)
```

---

## 5. Visualizations

### 5.1 The Bayesian Clinical Volcano Plot
```r
plot.pnp_individual_certainty(
  x,
  mid_select = NULL,
  point_alpha = 0.85,
  interval_alpha = 0.5,
  ...
)
```
1. **Axes & Geometry:**
   - **x-axis:** Individual posterior `median`.
   - **Horizontal Error Bars:** Representing `ci_low` to `ci_high`.
   - **y-axis:** `directed_certainty` mapped from `0%` to `100%` (using `scales::label_percent()`).
2. **ROPE Shading:**
   - Use `geom_rope()` to cleanly shade the vertical region $[-m, +m]$. 
3. **Automatic Directional Annotations:**
   - Read `x$meta$better_is` to map the directional arrows at the bottom/top of the plot:
     - `better_is == "lower"`: Left is `"\u2190 Favors Improvement"`, Right is `"Favors Deterioration \u2192"`.
     - `better_is == "higher"`: Left is `"\u2190 Favors Deterioration"`, Right is `"Favors Improvement \u2192"`.

### 5.2 The Sensitivity Grid Plot
Design a dedicated visualization function that accepts the summarized tibble output generated by the existing `get_clinical_sensitivity_grid()` function.

```r
plot_clinical_sensitivity(
  grid_data,
  color_palette = "viridis",
  ...
)
```
- **Inputs:** `grid_data` must contain `rope`, `threshold`, and `prop_certain` (as output by `get_clinical_sensitivity_grid()`).
- **Mapping:** 
  - `x = rope`
  - `y = prop_certain` (formatted as percentages).
  - `color = as.factor(threshold)` and `group = as.factor(threshold)`.
- **Geoms:** Use `geom_line()` and `geom_point()` to show how the proportion of the cohort meeting clinical significance (`prop_certain`) decays as the Minimal Important Difference (`rope`) widens, across varying required certainty `threshold` levels (e.g., 80%, 90%, 95%).
- **Theme:** Minimal, publication-ready styling ensuring clear differentiation between the threshold lines.
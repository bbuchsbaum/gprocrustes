# RFE-0001: First-class entity correspondence and column-space contracts in `multidesign`

| Field | Value |
|---|---|
| Status | Draft |
| Date | 2026-08-20 |
| Target package | [`multidesign`](https://github.com/bbuchsbaum/multidesign) (`~/code/multidesign`) |
| Motivating consumers | [`gprocrustes`](https://github.com/bbuchsbaum/gprocrustes) (this repo), [`manifoldalign`](https://github.com/bbuchsbaum/manifoldalign) |
| Related | `multivarious` projector generics; `gprocrustes` `proc_data` / `row_map` / cell masks |
| Author | Assessment from the gprocrustes 1.0 work |

This is a request for enhancement against **`multidesign`**, filed here because the need was discovered while deciding whether `gprocrustes` can reuse `hyperdesign` as input. Implementation belongs in `multidesign`. `gprocrustes` remains an adapter client.

---

## 1. Summary

Add an **opt-in contract** to `hyperdesign` (and a small slot on `multidesign`) that names:

1. a **global entity key** — one design column whose *values* identify the same row-thing across blocks;
2. a **column-space role** — whether columns are a shared Euclidean space \(\mathbb{R}^d\) or block-specific features;
3. later, an optional **cell observation mask**.

Do **not** add a flipped storage mode (`d × n`, features-by-observations). Rows stay the design axis. Left-action solvers continue to transpose at their own boundary.

Today `hyperdesign` only records `common_vars`: the intersection of design *column names*. That is metadata schema, not correspondence. `manifoldalign::generalized_procrustes.hyperdesign` already reconstructs correspondence by hand from a `y` / `task` column and then transposes. This RFE asks `multidesign` to own the first half of that, so every consumer does not re-derive it.

---

## 2. Problem

### 2.1 What `multidesign` is good at

A `multidesign` is a matrix plus a row design tibble plus optional `column_design`:

- `x`: rows = observations, columns = variables
- `design`: one row per observation, including a reserved `.index`
- `column_design`: one row per variable

A `hyperdesign` is a named list of those blocks (subjects, sessions, modalities) with:

- `attr(, "hdes")` — per-block `nr`, `nxvar`, row/col start–end
- `attr(, "common_vars")` — `Reduce(intersect, names(design))` minus `.index`

All verbs assume the **design axis is the rows**: `subset()`, `split()`, `summarize_by()`, `fold_over()`, `cv_rows()`. Column verbs (`select_variables()`, `column_design()`) cut the other axis. That API is correct and should not be inverted.

### 2.2 What three consumers need that it does not say

| Consumer | Rows mean | Columns mean | Cross-block join |
|---|---|---|---|
| Typical `multidesign` / KEMA / SSMA | trials, scans, time points (often i.i.d.-ish within a design cell) | block-specific features; `ncol` may differ | shared *factor names* (`condition`, `run`); alignment is learned |
| `manifoldalign` GPA | tasks / stimuli with a global label | features in a common \(d\); then **transposed** to \(d × n\) for \(O^\top A\) | a designated design column (`task`) |
| `gprocrustes` | corresponding entities (landmarks, vertices, stimuli) | dimensions of a common \(\mathbb{R}^d\); **right** action \(sXR+\mathbf{1}t\) | `ids` / `row_map` / `global_ids`; missing **rows** ≠ missing **cells** |

The storage layout of (entity × dimension) already matches (observation × variable). The mismatch is **semantics**, not `nrow`/`ncol`:

1. **Correspondence is not first-class.** Sharing a column named `task` is not the same as sharing the *value* `"landmark_17"` across blocks. `common_vars` cannot answer “which global entities does block 2 observe?” Partial overlap, the GPA default, is invisible.
2. **Column space is not declared.** `hyperdesign` allows unequal `ncol` (correct for multiblock embeddings) and equal `ncol` (required for GPA / `as_multidesign()` collapse). Nothing records which was intended. `as_multidesign()` *errors* if `ncol` differs, but that is a collapse constraint, not a space contract.
3. **Cell missingness does not exist.** Dropping a row can express “this landmark is absent.” It cannot express “this landmark is present and coordinate 2 is missing.” `gprocrustes` treats those as different mathematical objects. `NA` in `x` is not an acceptable substitute (implicit zeros / `NA` are not incidence).
4. **Collapse and correspondence are opposite operations.** `as_multidesign(hd)` *row-stacks* blocks. GPA *aligns* blocks onto a global entity index. Both are valid. Only stacking is implemented, and it is easy to reach for when the user meant join-by-id.

### 2.3 Why a storage-flip mode is the wrong fix

Left-action callers (`QX`, features × observations) are real. Putting that layout *inside* `multidesign` would:

- invert every row verb (`subset(md, task == "A")` becomes a column select);
- force `design` and `column_design` to swap or lie;
- break `manifoldalign`, KEMA, SSMA, and `cv_alignment_rows`, which all assume rows carry design;
- hide the orientation change that `gprocrustes` requires to be explicit at an adapter boundary.

`multiblock` already has `cstacked` / `rstacked`. That is about how matrices are bound, not which side a group acts on. Do not reuse it as an action flag.

**Action convention stays off the data object.**

---

## 3. Goals

**G1.** A `hyperdesign` may declare one design column as the **global entity key**. From that column, the package can build a correspondence object: global id universe, per-block `ids`, per-block integer `row_map`, and an incidence / overlap description.

**G2.** A `hyperdesign` may declare a **column-space role**: `common` (all blocks live in the same \(\mathbb{R}^d\); `ncol` must agree; `column_design` should be compatible) versus `block` (columns are per-block features; `ncol` may differ). Unset remains legal and means “unspecified,” i.e. today’s behaviour.

**G3.** Existing objects, verbs, and callers keep working with both attributes unset. No required arguments on `hyperdesign()` / `multidesign()`.

**G4.** Correspondence is derived from **values**, not from positional row index, unless the caller explicitly asks for positional correspondence (equal `nrow`, no id column). The assumption used must be recorded (`explicit_ids` vs `positional`), matching `gprocrustes::proc_data()`.

**G5.** Duplicate entity keys inside one block are an error unless an aggregation rule is supplied. This is already a `gprocrustes` law and should be the same law here.

**G6.** Optional later: a `cells` mask on `multidesign`, aligned to `x`, logical, no `NA`. Missing cells are not missing rows. `NA` in `x` is never interpreted as missingness by `multidesign` itself.

**G7.** Document that `fold_over(hd, <id>)` holds out **entities** when rows are entities. That is useful for landmark-holdout CV. It is not a new CV engine and must not be confused with trial-wise `fold_over(hd, run)`.

---

## 4. Non-goals

- Do not add `orientation = "variables_by_observations"` or store \(X^\top\) as a first-class mode.
- Do not put left vs right action, polar factors, gauges, or optimality certificates on the data class.
- Do not make `multidesign` depend on `gprocrustes`.
- Do not implement GPA, GPM, Gower, IRLS, or LBW.
- Do not turn `hyperdesign` into a `multivarious` projector, and do not implement `project` / `scores` / `ncomp` / `sdev` on it.
- Do not treat `NA` in `x` as a missingness model.
- Do not change `common_vars` to mean entity ids.
- Do not change `as_multidesign()` default from row-stack to id-align (see §8; a *new* function may align).
- Do not require unique global coverage (every entity in every block). Partial overlap is the point.
- Do not solve unknown correspondence (ICP, OT, Hungarian on features). Known identity is the contract, as in `gprocrustes`.
- Do not replace `gprocrustes::proc_data`. The engine-facing object stays there.

---

## 5. Current behaviour (baseline)

Relevant constructors and attributes, as of the local checkout `~/code/multidesign`:

```r
multidesign.matrix <- function(x, y, column_design = NULL, ...)
# nrow(x) == nrow(y); ncol(x) == nrow(column_design)
# design$.index <- seq_len(nrow(x))

hyperdesign.list <- function(x, block_names = NULL)
# all elements inherit "multidesign"
# attr(, "hdes")$nr, $nxvar, row_start/end, col_start/end
# attr(, "common_vars") <- intersect(names(design)) minus .index
```

`as_multidesign.hyperdesign` requires equal `ncol` and identical `column_design`, then `rbind`s `x` and `bind_rows` the common design columns. That is **stack subjects**, not **match landmarks**.

`manifoldalign` (`R/genprocrustes.R`, `generalized_procrustes.hyperdesign`):

1. optional `multivarious` preprocess (`center()`, `init_transform`);
2. `design[[y]]` as task labels;
3. `t(x)` so the solver sees features × tasks;
4. integer map over `unique(unlist(labels))`.

`gprocrustes` (`R/data.R`, `proc_data()`):

- `views`, `ids`, `observed` (row mask), `cells` (entity × dimension mask);
- `.gproc_correspondence()`: all views have ids (or rownames) → `explicit_ids`; none do and `nrow` equal → `positional`; mixed → error;
- duplicate ids in a view → `duplicate_entity_ids`.

The RFE is: move the *correspondence construction* (step 2 of manifoldalign; the id half of `proc_data`) into `multidesign`, and declare column space, so adapters become mechanical.

---

## 6. Proposed contract

### 6.1 Attributes on `hyperdesign`

| Attribute | Type | Default | Meaning |
|---|---|---|---|
| `common_vars` | character | as now | shared design *column names* (unchanged) |
| `id` | `NULL` or length-1 character | `NULL` | name of the design column that is the global entity key |
| `space` | `NULL`, `"common"`, or `"block"` | `NULL` | column-space role |
| `correspondence_assumption` | `NULL`, `"explicit_ids"`, or `"positional"` | `NULL` | how row identity was resolved; set when `id` is set or when `positional = TRUE` |

Suggested constructor extension (all new arguments optional):

```r
hyperdesign(
  x,
  block_names = NULL,
  id = NULL,
  space = NULL,
  positional = FALSE
)
```

Validation at construct time:

- If `id` is not `NULL`:
  - `id` must be in every block’s `design` (hence also in `common_vars`).
  - `id` must not be `.index` or `.orig_index`.
  - Within each block, `design[[id]]` has length `nrow(x)`, no `NA`, and no duplicates unless `aggregate` is later supplied (phase 2).
  - Coerce keys to character for matching (`"17"` and `17` join). Preserve the original column type on the design tibble.
- If `positional = TRUE`:
  - `id` must be `NULL`.
  - All blocks have equal `nrow`.
  - Record `correspondence_assumption = "positional"`.
- If `space = "common"`:
  - all `ncol(x)` equal;
  - `column_design` identical or *compatible* (same `nrow`, same `.index` if present; extra metadata columns may be a later relaxation — start strict: `identical()` as `as_multidesign()` already does).
- If `space = "block"`:
  - unequal `ncol` is allowed;
  - `as_multidesign()` and any “shared \(\mathbb{R}^d\)” helper must refuse.
- If `space` is `NULL`: no new checks (today’s behaviour).

`hyperdesign()` of a list that already carries `id` / `space` on the blocks: the list-level attributes win if supplied; otherwise copy a unanimous block-level setting; if blocks disagree, error.

### 6.2 Optional slot on `multidesign`

```r
multidesign(x, y, column_design = NULL, cells = NULL, ...)
```

- `cells`: `NULL` or a logical matrix with `dim(cells) == dim(x)`, no `NA`.
- A row whose `cells` row is all `FALSE` is still a row of the design (the entity is listed but fully unobserved in coordinates). That is not the same as omitting the row. Document this; do not auto-drop.
- Verbs that subset rows (`subset`, `fold_over`, `cv_rows`) subset `cells` on the same indices. Verbs that select columns (`select_variables`) subset `cells` columns.
- `summarize_by()` in phase 1 errors if `cells` is present (aggregation of partially observed cells is a consumer problem). Phase 2 may accept an explicit `na.rm`-style rule that still does not treat `NA` in `x` as the mask.

Do not add `id` as a required column on a single `multidesign`. A lone block has no cross-block join. The key lives on the `hyperdesign`.

### 6.3 Correspondence object

New S3 class `md_correspondence` (name bikeshed, see §12), returned by `correspondence(hd)`:

```r
correspondence(x, ...)  # generic

# hyperdesign method
# requires id or positional assumption; otherwise errors with a pointer
# to hyperdesign(..., id = )

list(
  id             = "landmark",          # or NULL if positional
  assumption     = "explicit_ids",      # or "positional"
  global_ids     = c("a", "b", "c"),    # character, union, first-seen order
  ids            = list(                # per block, character, length nrow
    specimen1 = c("a", "c"),
    specimen2 = c("a", "b", "c")
  ),
  row_map        = list(                # per block, integer index into global_ids
    specimen1 = c(1L, 3L),
    specimen2 = c(1L, 2L, 3L)
  ),
  n_global       = 3L,
  n_observed     = c(specimen1 = 2L, specimen2 = 3L),
  overlap        = <see below>
)
class = c("md_correspondence", "list")
```

**Overlap.** Minimum useful form, enough for a disconnected-component check without pulling in `gprocrustes`:

- `incidence`: \(n_{\mathrm{global}} \times K\) logical (or sparse `dgCMatrix` 0/1) — entity \(j\) observed in block \(k\);
- `pair_n`: \(K \times K\) integer — number of shared global ids between blocks;
- `connected`: logical — the undirected graph with an edge when `pair_n[i,j] > 0` is connected.

Do not claim this is the Ling overlap graph used for GPM certificates. It is the design-level incidence graph. Consumers may build a richer graph from it.

**Order of `global_ids`.** First-seen: walk blocks in list order, append unseen keys. Document this. Do not sort unless the user passes `sort_ids = TRUE`. Sorting would silently permute consensus rows for every caller.

### 6.4 Helpers

```r
has_correspondence(hd)           # id set or positional assumption recorded
entity_id(hd)                    # the id name, or NULL
column_space(hd)                 # "common", "block", or NULL

align_by_id(hd, fill = NA_real_) # NEW; see §8. Not the default of as_multidesign()
```

`xdata()`, `design()`, `column_design()` unchanged.

Print methods should show `id`, `space`, \(n_{\mathrm{global}}\), and whether the overlap graph is connected, when those attributes are set.

---

## 7. Worked examples

### 7.1 Partial landmark overlap (GPA input)

```r
library(multidesign)

# Rows = landmarks, columns = coordinates. Same layout as today's multidesign.
s1 <- multidesign(
  matrix(rnorm(10 * 3), 10, 3),
  data.frame(landmark = paste0("L", 1:10), curve = rep(c("jaw", "orbit"), 5)),
  data.frame(axis = c("x", "y", "z"))
)
s2 <- multidesign(
  matrix(rnorm(8 * 3), 8, 3),
  data.frame(landmark = paste0("L", c(1:6, 11, 12)), curve = rep(c("jaw", "orbit"), 4)),
  data.frame(axis = c("x", "y", "z"))
)

hd <- hyperdesign(
  list(specimen1 = s1, specimen2 = s2),
  id = "landmark",
  space = "common"
)

cr <- correspondence(hd)
cr$global_ids          # L1..L12 (first-seen)
cr$row_map$specimen2   # indices of L1..L6, L11, L12 in that universe
cr$overlap$connected   # TRUE if the two views share at least one landmark
```

`gprocrustes` adapter (not part of this RFE’s implementation, shown as the acceptance consumer):

```r
as_proc_data.hyperdesign <- function(x, ...) {
  cr <- correspondence(x)
  views <- lapply(x, function(d) d$x)
  cells <- lapply(x, function(d) d$cells)  # may be all NULL
  proc_data(views, ids = cr$ids, cells = cells)
}
```

No transpose. If a left-action solver wants \(d × n\), it transposes after this, visibly.

### 7.2 Task-labelled domains (today’s manifoldalign pattern)

```r
d1 <- multidesign(
  matrix(rnorm(5 * 10), 5, 10),
  data.frame(task = factor(c("A", "B", "C", "D", "E")))
)
d2 <- multidesign(
  matrix(rnorm(4 * 10), 4, 10),
  data.frame(task = factor(c("A", "C", "D", "F")))
)
hd <- hyperdesign(list(domain1 = d1, domain2 = d2), id = "task", space = "common")

# manifoldalign GPA would become:
#   labels <- correspondence(hd)$ids
#   A_list <- lapply(xdata(hd), t)    # still their flip, still their solver
```

The `y` argument on `generalized_procrustes(hd, task)` can default to `entity_id(hd)` when set.

### 7.3 Unchanged KEMA-style hyperdesign

```r
hd <- hyperdesign(list(subj1 = md1, subj2 = md2))
# id = NULL, space = NULL
# correspondence(hd) errors
# fold_over(hd, run) unchanged
```

If someone later sets `space = "block"` to be honest about unequal feature counts, that is documentation plus a refusal from `align_by_id()` / `as_multidesign()`.

### 7.4 Positional correspondence

```r
hd <- hyperdesign(list(md1, md2, md3), positional = TRUE, space = "common")
correspondence(hd)$assumption  # "positional"
# global_ids = as.character(seq_len(nrow))
```

Refuse `positional = TRUE` with unequal `nrow`, even if a human “knows” they line up after dropping rows.

---

## 8. `as_multidesign()` versus `align_by_id()`

These must not be unified.

| Function | Geometry | Requires | Result `nrow` |
|---|---|---|---|
| `as_multidesign(hd)` | row-stack blocks | equal `ncol`, compatible `column_design` | \(\sum_k n_k\) |
| `align_by_id(hd)` | place each block’s rows onto the global entity index | `id` or positional; `space = "common"` | \(n_{\mathrm{global}}\) per block, then typically a list or a 3-way array, **not** a single stacked `multidesign` |

`align_by_id()` proposed return (phase 2 is fine):

```r
# list of matrices, each n_global × d, missing entities filled with `fill`
# plus the correspondence object
# OR an n_global × d × K array when the caller asks
```

Filling with `0` is forbidden as a default (`gprocrustes`: implicit zeros are never missing). Default `fill = NA_real_` is only a *presentation* for human inspection. Any solver adapter must use the incidence map, not the filled array, as the missingness model.

Phase 1 can ship `correspondence()` only and leave `align_by_id()` to consumers. Phase 2 ships the helper so people stop abusing `as_multidesign()` for GPA.

---

## 9. Interaction with existing verbs

| Verb | With `id` / `space` unset | With contract set |
|---|---|---|
| `subset(hd, expr)` | as now | rebuild `hyperdesign()` so `id` / `space` are preserved; correspondence is recomputed (global universe may shrink) |
| `select_variables(hd, expr)` | as now | preserve `id`; if `space = "common"`, all blocks stay equal-`ncol` or error |
| `fold_over(hd, run)` | as now | as now; `run` is not the entity key |
| `fold_over(hd, landmark)` | works today as “hold out rows with that factor level” | **document** that this holds out entities; still not Bai–Bartoli math; `gprocrustes::cross_validate()` stays the LBW engine |
| `cv_rows()` | as now | row indices remain local; `preserve_row_ids` already exists — if `id` is set, also preserve `design[[id]]` (it already would, as an ordinary column) |
| `summarize_by(hd, condition)` | `colMeans` per design cell | if `id` is set, error unless the user is clearly aggregating *replicates of the same entity* via an explicit `aggregate` rule (phase 2). Silent `colMeans` over landmarks grouped by `curve` is a different analysis and remains allowed when grouping vars ≠ `id` |
| `init_transform(hd, preproc)` | per-block preprocess | allowed; does not change `id`. If preprocess changes `ncol`, `space = "common"` must be revalidated |
| `reduce.multidesign` | PCA on columns | a reduced block is `space = "block"` unless the user reasserts `common` after reducing every block to the same \(d\) |
| `block_indices(..., byrow =)` | as now | unchanged; still about concatenated index ranges, not correspondence |

Reconstruction after subset: `hyperdesign(filtered_blocks, id = entity_id(old), space = column_space(old))`, not a hand-copied attribute list that can drift.

---

## 10. Compatibility and dependencies

- **No new hard dependencies.** Correspondence can be built with base R. Sparse incidence may use `Matrix` if `multidesign` already imports it; otherwise a dense logical matrix is enough for phase 1 (typical \(n_{\mathrm{global}}\) is landmarks or tasks, not voxels).
- **Do not import `gprocrustes` or `manifoldalign`.**
- `multivarious` remains what it is today (`reduce()`, `init_transform()`). This RFE does not expand that surface.
- S3 class names and `common_vars` stay. Old `hyperdesign` objects without the new attributes remain valid.
- If attributes are stored with `attr(hd, "id")`, avoid clashing with `as_multidesign(hd, .id = )`, which today means “name of the column that *labels the source block* after stacking.” Different concept. Prefer `attr(hd, "entity_id")` internally if `"id"` is too overloaded. Public argument can still be `id =`.

---

## 11. Suggested API (phase 1 vs phase 2)

### Phase 1 — contract and query (this RFE’s minimum)

- `hyperdesign(..., id = NULL, space = NULL, positional = FALSE)`
- `correspondence.hyperdesign`
- `has_correspondence()`, `entity_id()`, `column_space()`
- preserve attributes through `subset.hyperdesign` and `select_variables.hyperdesign`
- `print.hyperdesign` shows the contract when set
- tests in §13
- documentation: a short section on “rows as entities” vs “rows as trials”; explicit warning not to transpose the stored `x`

### Phase 2 — masks and align helper

- `multidesign(..., cells = NULL)` and row/column verb propagation
- `align_by_id()`
- optional `aggregate` for duplicate keys (`mean`, or a user function) — default remains error
- optional `sort_ids`

### Phase 3 — out of scope unless a consumer asks again

- `gpa.hyperdesign` / `as_proc_data.hyperdesign` in **`gprocrustes`**, Suggests: `multidesign`
- `generalized_procrustes` reading `entity_id(hd)` in **`manifoldalign`**
- sparse incidence, disconnected-component diagnostics beyond `overlap$connected`

Phase 3 is listed so `multidesign` does not grow solvers “to be helpful.”

---

## 12. Alternatives considered

**A. Do nothing; each consumer keeps parsing `design[[y]]`.**
Rejected as the long-term state. `manifoldalign` and `gprocrustes` would duplicate id validation, duplicate-key errors, and overlap. They already disagree on action convention; they should not also disagree on what an id is.

**B. Flipped `multidesign` / `orientation =` storage flag.**
Rejected. See §2.3. The tidy API is row-centric. Left-action is a solver boundary.

**C. New class `entitydesign` / `shapedesign`.**
Honest, but splits the ecosystem. GPA data *is* a `multidesign` whose row key is an entity and whose columns are a shared \(\mathbb{R}^d`. A flag-plus-validation on `hyperdesign` is enough. Revisit only if `cells` + correspondence force every method to grow an unreadable branch.

**D. Encode correspondence as rownames only.**
Insufficient. rownames are easy to lose under `subset` / `as.matrix`; they do not travel through `design`; mixed “some blocks have rownames” is already an error in `proc_data`. A design column is the right place (and is what `manifoldalign` uses).

**E. Let `common_vars` grow a special name like `.entity`.**
Too magical. `common_vars` means “these columns exist in every design.” An entity key is one of those columns plus a join interpretation. Keep the attributes separate.

**F. Make `as_multidesign()` align by id when `id` is set.**
Rejected. Silent change of geometry (stack vs join) will corrupt existing stacking workflows the first time someone adds `id =` for documentation. New name: `align_by_id()`.

**G. Put `id` on each `multidesign` and intersect at `hyperdesign()` time.**
Acceptable variant. Slightly more verbose (`multidesign(..., id = "landmark")` on every block) and then `hyperdesign()` checks they match. Prefer list-level `id` for less repetition; allow block-level as a copy source (see §6.1).

---

## 13. Acceptance tests

All of these are `multidesign` tests. None require loading `gprocrustes`.

1. **Unset contract is today’s object.** `hyperdesign(list(d1, d2))` has `entity_id() == NULL`, `column_space() == NULL`, and existing tests still pass.
2. **`id` must exist on every block.** Missing column → construction error, not a later `correspondence()` surprise.
3. **Duplicate keys in one block → error** (phase 1).
4. **`NA` in the id column → error.**
5. **Partial overlap.** Two blocks, ids `{A,B,C}` and `{B,C,D}` → `global_ids` length 4, `row_map` correct, `pair_n[1,2] == 2`, `connected == TRUE`.
6. **Disconnected overlap.** Blocks `{A,B}` and `{C,D}` → `connected == FALSE`. Construction succeeds; the flag is data, not an error. (Consumers decide whether that is `invalid_problem`.)
7. **First-seen order.** Block1 `{B,A}`, block2 `{A,C}` → `global_ids == c("B","A","C")`.
8. **Type coercion.** Integer `1:3` and character `c("1","2")` join; no leftover unmatched `"1"` vs `1`.
9. **`space = "common"` rejects unequal `ncol`.**
10. **`space = "block"` allows unequal `ncol` and makes `as_multidesign()` refuse with a message that names `space`.**
11. **`positional = TRUE` requires equal `nrow` and `id = NULL`.**
12. **`subset` preserves `id` and `space` and recomputes correspondence.** Subsetting away all shared ids can yield `connected == FALSE`.
13. **`select_variables` on `space = "common"` keeps equal `ncol` or errors.**
14. **`fold_over` / `cv_rows` still run** on a contracted `hyperdesign`; assessment rows keep the id column.
15. **`cells` (phase 2):** wrong `dim` errors; `NA` in mask errors; `subset` on rows subsets `cells`; a row of all-`FALSE` is retained.
16. **No `NA` in `x` is treated as a mask** by `correspondence()` or `align_by_id()`.
17. **`align_by_id` (phase 2)** default fill is `NA`, never `0`; output rows follow `global_ids`.

---

## 14. Documentation requirements

When this lands in `multidesign`:

- README / vignette: one paragraph that rows may be trials *or* corresponding entities, and that this is a **role of the row key**, not a second matrix orientation.
- Explicit sentence: **do not store features × observations in `x` to please a left-action formula.** Transpose in the solver.
- Contrast `common_vars` (schema) with `id` (join).
- Contrast `as_multidesign()` (stack) with `align_by_id()` (join).
- Contrast `space = "common"` (GPA, shared \(\mathbb{R}^d\)) with `space = "block"` (multiblock embeddings).
- Point at `gprocrustes` for transforms, gauges, certificates, and cell-aware *fitting*; this package only carries the mask.
- `fold_over(hd, <entity id>)` is “hold out those rows,” not “the LBW smoothness CV of Bai–Bartoli.”

---

## 15. Implications for `gprocrustes` (consumer, not this RFE’s patch)

Once phase 1 exists, Milestone 5 adapters stay small:

```r
# Suggests: multidesign
gpa.hyperdesign <- function(data, ..., id = entity_id(data)) {
  if (is.null(id) && !has_correspondence(data)) {
    stop("hyperdesign needs id= or positional correspondence")
  }
  if (!identical(column_space(data), "common")) {
    stop("gpa() requires space = \"common\"")
  }
  gpa(as_proc_data(data), ...)
}
```

Still out of scope for `multidesign`: typed transforms, solver dispatch, certificates, \(\Sigma_S\) vs \(\Sigma_M\), broom methods, lazy aligned stores.

`gprocrustes` will **not** take `multidesign` as a hard Import for 1.0. The freeze already names `proc_data()` / `proc_from_array()` as the data constructors.

---

## 16. Implications for `manifoldalign` (consumer)

`generalized_procrustes.hyperdesign(data, y, ...)` can:

- if `missing(y)` and `entity_id(data)` is set, use that column;
- take `correspondence(data)$ids` instead of a local `purrr::map` + `unique(unlist(...))`;
- keep `t(x)` and GPM **inside manifoldalign** (left action, their object).

That is a later cleanup, not a `multidesign` blocker.

---

## 17. Open questions (need a human decision in `multidesign`)

1. **Public attribute name:** `id` vs `entity_id` vs `key`. `id` is short and matches `gprocrustes`; it collides with `as_multidesign(.id=)` in prose.
2. **`space` vocabulary:** `"common"` / `"block"` vs `"dimension"` / `"feature"` vs `"shared"` / `"separate"`. Prefer `"common"` / `"block"` — it matches how people already talk about hyperdesign blocks.
3. **Strict `column_design` identity** for `space = "common"`, or merely equal `ncol`? Start strict (today’s `as_multidesign` rule). Relax later if someone has per-block column annotations on the same axes.
4. **Phase 2 `cells` on `multidesign` vs only on `hyperdesign`.** Prefer the single-block slot so subset methods have one place to update.
5. **Should construction fail on `connected == FALSE`?** No. Report it. Disconnected GPA is a solver / compiler decision (`gpa(..., allow_disconnected = )`).
6. **Integer vs character keys in `row_map`.** Integer index into `global_ids` (as `proc_data`) is the right export. Do not export a list of factors that drop unused levels differently per block.

---

## 18. Suggested issue text for the `multidesign` tracker

> **Feature:** Opt-in entity correspondence and column-space contract on `hyperdesign`.
>
> `common_vars` records shared design *column names*. GPA and other known-correspondence methods need shared *entity values* and a declaration that columns are a common \(\mathbb{R}^d\). Add optional `id`, `space`, and `correspondence()`. Do not add a transposed storage mode. Do not change `as_multidesign()` from row-stack to join.
>
> Details: gprocrustes `docs/rfe/0001-multidesign-correspondence-and-space.md`.
>
> Acceptance: phase 1 tests 1–14 in that RFE.

---

## 19. Decision requested

From `multidesign` maintainers:

1. Accept phase 1 as specified (§11), or
2. Accept phase 1 with a different naming choice from §17, or
3. Reject in favour of a new class (alternative C), or
4. Reject — consumers keep parsing design columns themselves.

Until that lands, `gprocrustes` keeps `proc_data` as the only solver-facing type and does not depend on `multidesign`.

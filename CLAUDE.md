# LCDmix — working context

Read this before touching anything. It explains what is being built, what must not break, and
which mistakes have already been made here.

## The project

LCDmix is a mixture-of-experts regression model: covariate-dependent multinomial-logistic gating
plus **log-concave** (nonparametric) expert error densities, fitted by an EM-type algorithm with an
ℓ₁ penalty on both gating and expert coefficients. It is implemented as an R package built with
`litr` — **the package is generated from `create-LCDmix.Rmd`**, so that Rmd is the source of truth.
Editing files under `LCDmix/R/` directly is wrong; they are build output.

A paper is being prepared for *Computational Statistics & Data Analysis*. Two planning documents in
the sibling folder `../` carry the full reasoning:

- `LCDmix_plan_B_sprint.md` — the 10-day plan to submission (self-contained)
- `LCDmix_ablation_design.md` — the experiment this code work exists to support
- `n_outside_diagnosis.md` — why the out-of-support counter matters and what it measures

Read `LCDmix_ablation_design.md` before proposing changes to the M-step. It explains the
statistical reasoning, and changes that look like cleanups can silently invalidate the experiment.

## What is being built right now

The paper claims, four times, that the feasibility-aware LP update for the expert coefficients
"improves numerical stability relative to unconstrained updates." **No experiment has ever tested
that.** The work in progress is a minimal ablation: the same EM, run twice, differing only in the
θ-update — the LP versus an unconstrained quasi-Newton solve.

Key statistical fact that motivates the design: `ĝ_k` is piecewise linear and concave, and composing
it with an affine map in θ keeps it so. **The θ-block objective is therefore a concave
piecewise-linear program.** The LP solves it exactly; a quasi-Newton method is structurally
ill-suited, because the optimum sits at a vertex where the gradient is undefined. The experiment
measures how much that costs, it does not establish the point.

## Current state

Implemented and working:

- `mstep_theta()` dispatches on a new `update = c("lp", "optim")` argument
- `mstep_theta_optim()` — the comparison arm; L-BFGS-B over `(θ⁺, θ⁻) ≥ 0`, same objective and same
  `N_total * lambda_theta` penalty scaling as `mstep_theta_lp()`
- `make_logdens_ext()` — evaluates `ĝ_k` with **linear extension beyond the support**, using the
  boundary slopes. This is deliberate: with the honest `-Inf` convention `optim` stalls on its first
  step outside the support and the comparison proves nothing. Do not "fix" this to return `-Inf`.
- `comp_Q()` gained appended optional args `Y_bin = NULL, intercepts = NULL`; when both are supplied
  it recomputes residuals from `(intercepts, slopes)` rather than trusting the `residuals` argument
- `comp_Q()` returns a **scalar** carrying diagnostics as attributes: `n_eval`, `n_outside`,
  `mass_outside`, `max_over`
- `iteration()` records `Q_every`, `n_outside_every`, `lp_check_every`, `theta_diag`, and propagates
  `error` / `failed_iter` from its `tryCatch`

Verified:

- The LP keeps residuals inside `[L_k, U_k]` on the bins it constrained. `lp_check_every` reports at
  most one bin with `max_over ≈ 4.44e-16` (= `2 * .Machine$double.eps`) — floating-point noise from
  `θ^(m)` sitting exactly on the boundary, since `L = min(res_k)`. This is the premise the whole
  ablation rests on and it holds.

## Open bug — start here

On a run with `calc_Q_every = TRUE` that converged at iteration 6:

```
length(fit$iter$Q_every)      # 25  -> 1 + 6*4
fit$iter$iter_num             # 6
fit$iter$n_outside_every      # 24 entries -> 6*4
# 0 0 1 0 | 0 0 0 0 | 0 0 3 0 | 0 0 2 0 | 0 0 3 0 | 0 0 1 0
```

There should be **five** gated/unconditional entries per iteration: after the E-step, after the α
update, after the θ update, after centering, and after the density update. Only four are recording.

Working hypothesis: the `if (calc_Q_every)` block at the **θ-update site** was dropped or altered
when the `lp_check` code was inserted. Verify by reading `iteration()` and counting the
`Q_every <- c(Q_every, Q_new)` sites; there must be exactly four inside `if (calc_Q_every)` plus one
unconditional after the density update.

This matters for interpretation, not just tidiness. If the hypothesis is right, the variable entries
above are the **centering step**, not the θ update — and "centering pushes residuals off the
support" is a substantive finding about the algorithm, so the slot labelling has to be certain
before it is written up.

Add a cheap assertion so this cannot drift again:

```r
if (calc_Q_every) stopifnot(length(Q_every) == 1L + 5L * i)
```

## Hard invariants

Breaking any of these invalidates saved results or the experiment.

1. **Do not modify the `comp_Q()` call after the density update.** It produces the `Q` that drives
   the stopping rule and the revert-and-break path. It is deliberately left with positional
   arguments and no recomputation so the convergence sequence stays bit-identical to the fits stored
   in `refit_bests.rds` and `cv_summ.rds`. Those saved objects are load-bearing for the paper's
   analysis and must remain reproducible.

2. **New function arguments are appended, never inserted.** `main()` calls `iteration()` positionally
   with twelve arguments; `comp_Q()` is called positionally at five sites. Inserting a parameter
   mid-signature does not error — it silently misbinds and produces wrong numbers.

3. **`mstep_theta_lp()` is called with `idx = idx_old`, deliberately.** That keeps `L, U` consistent
   with the `ĝ_k` those residuals were fitted to, which is what makes the feasible set non-empty at
   the current iterate. Do not "fix" it to `idx_new`. (The manuscript currently says `idx_new`; the
   manuscript is what gets corrected, not the code.)

4. **`comp_Q()` returns a scalar.** Diagnostics are attributes. Returning a list breaks every call
   site, since `Q_new` is consumed by `c(Q, Q_new)`.

5. **Do not change the algorithm mid-sprint.** Folding the centering constraint into the LP is a good
   idea and is explicitly deferred to the revision round — see the deferred register in
   `LCDmix_plan_B_sprint.md`.

## Mistakes already made here — do not repeat

Every bug in this codebase so far has been a wiring error, not a logic error. The pattern is a
snippet applied without its declaration:

| Symptom | Actual cause |
|---|---|
| `object 'resi_new' not found` | `update = update` used in `iteration()` before `update` was a parameter; R resolved it to `stats::update`, `match.arg()` threw, `tryCatch` swallowed it, and the return statement outside the `tryCatch` then referenced an unassigned variable |
| `object 'max_over' not found` | counter used in `comp_Q()` without `max_over <- 0` beside the other initializers |
| `argument 7 is empty` | trailing comma in `main()`'s partial-fit `return(list(...))` — pre-existing dead code, made reachable for the first time by the new error propagation |
| `lp_check_every` empty | assigned in the loop, never declared or returned |

**When adding a counter or accumulator, add the declaration, the assignment, the `last_state` entry,
and the return entry in the same edit.** Grep for the name afterwards and confirm four hits.

Also note: `stopifnot(is.null(fit$iter$error))` passes vacuously if the `error` field is not in the
return list. `$` on a missing name gives `NULL`. Prefer `stopifnot("error" %in% names(fit$iter))`
alongside it.

## Verification

`calc_Q_every` defaults to **FALSE** in both `main()` and `iteration()`. Without it, only the
unconditional post-density entry records — and that one is trivially zero, because `g_new` is refit
to `resi_new`, so its support covers those residuals by construction. A vector of zeros from a
`calc_Q_every = FALSE` run means nothing.

```r
fit <- main(Y = Y_bin, X = X, biomass = bin_mass, K = 2,
            binned = TRUE, calc_Q_every = TRUE, debug = TRUE)

stopifnot("error" %in% names(fit$iter), is.null(fit$iter$error))
stopifnot(length(fit$iter$Q_every) == 1L + 5L * fit$iter$iter_num)   # currently FAILS

do.call(rbind, fit$iter$lp_check_every)   # n_out 0, or max_over < 1e-8
```

Use a tolerance when counting out-of-support bins — `1e-8` — or boundary-touching residuals are
miscounted as violations.

Backward-compatibility checks that have already passed and should keep passing:

```r
# comp_Q with the new args omitted must equal the old value (attributes are new)
all.equal(as.numeric(comp_Q(X, g, resi, th, al, idx, w, la, lt)), Q_reference)

# supplying residuals directly vs recomputing must agree when they match
resi_chk <- comp_resi(Y_bin, X, theta0, theta)
all.equal(as.numeric(comp_Q(X, g, resi_chk, theta, al, idx, w, la, lt)),
          as.numeric(comp_Q(X, g, resi_chk, theta, al, idx, w, la, lt,
                            Y_bin = Y_bin, intercepts = theta0)))
```

## Next steps

1. Fix the missing θ-site recording; add the length assertion.
2. Re-run and confirm five entries per iteration, then read slots 3 and 4 — the θ update and the
   centering step — with `max_over` alongside the counts to separate real overshoot from
   floating-point noise.
3. Smoke test: one replication, both arms, and confirm the **LP arm reproduces the existing fit
   exactly** after the refactor. If it does not, the threading is wrong and every later number is
   suspect. Do not skip this.
4. Thread `update` through `cv_lcd_onejob()` and `refit_onejob()`.
5. Launch the ablation: one configuration (s_α = 10, θ₂₀ = 5), 20 replications, both arms, paired on
   seeds, penalties fixed at the CV-selected LCDmix values. Log per run: convergence/failure,
   objective decreases, iterations, wall-clock, final penalized log-likelihood, and the
   out-of-support counts.

Metrics and the reporting table are specified in `../LCDmix_ablation_design.md` §5.

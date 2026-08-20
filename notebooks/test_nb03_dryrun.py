"""Headless dry run for Notebook 03 (Bayesian optimization controller).

Exercises the core, GUI-free logic of cells 3, 4, 6, 7, 8, and 10 against
mock data and verifies the 10 checks listed in the build spec.
"""
import itertools
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid", context="paper", font_scale=1.3)

PLOTS_DIR = Path("/tmp/nb03_test_plots")
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

results = {}


def check(n, name, cond, actual=None, expected=None):
    status = "PASS" if cond else "FAIL"
    results[n] = status
    extra = ""
    if actual is not None or expected is not None:
        extra = f"  (actual={actual!r}, expected={expected!r})"
    print(f"CHECK {n}: {name} -> {status}{extra}")


# ---------------------------------------------------------------------------
# Mock configuration (mirrors Cell 2 defaults; no widgets/popups here)
# ---------------------------------------------------------------------------
candidate_pool = [55, 117, 119, 248, 249, 281, 308, 311, 312, 337, 361]
contact_scores = {
    55: 0.18, 117: 0.82, 119: 0.75, 248: 1.00, 249: 0.82,
    281: 0.50, 308: 0.78, 311: 0.82, 312: 0.85, 337: 0.84, 361: 0.60,
}
FORCE_INCLUDE_RESIDUES = [119, 312]
MIN_RESIDUES = 3
MAX_RESIDUES = 8
BATCH_SIZE = 5
XI = 0.01
USE_CONTACT_PRIOR = True
ACQUISITION_FUNC = "Expected Improvement"
CONVERGENCE_PATIENCE = 10
CONVERGENCE_THRESHOLD = 0.1

print("=" * 70)
print("NOTEBOOK 03 DRY RUN — headless core-logic verification")
print("=" * 70)

# ---------------------------------------------------------------------------
# Cell 3 logic: build combinatorial space
# ---------------------------------------------------------------------------
print("\n--- Cell 3: candidate pool & combinatorial space ---")

ALL_COMBINATIONS = []
for k in range(MIN_RESIDUES, MAX_RESIDUES + 1):
    for combo in itertools.combinations(candidate_pool, k):
        if all(r in combo for r in FORCE_INCLUDE_RESIDUES):
            ALL_COMBINATIONS.append(combo)


def encode(combo, pool=candidate_pool):
    combo_set = set(combo)
    return np.array([1.0 if r in combo_set else 0.0 for r in pool])


ALL_COMBINATIONS_ENCODED = np.array([encode(c) for c in ALL_COMBINATIONS])

raw_contact_sums = np.array([sum(contact_scores[r] for r in c) for c in ALL_COMBINATIONS])
max_contact_sum = raw_contact_sums.max() if len(raw_contact_sums) else 1.0
COMBO_CONTACT_SCORES = raw_contact_sums / max_contact_sum

print(f"Candidate pool: {candidate_pool}")
print(f"Force-include: {FORCE_INCLUDE_RESIDUES}")
print(f"Combination constraints: min={MIN_RESIDUES}, max={MAX_RESIDUES} residues")
print(f"Valid combinations: {len(ALL_COMBINATIONS)} (out of {2 ** len(candidate_pool)} total)")

fig, ax = plt.subplots(figsize=(7, 4))
cs_sorted = sorted(contact_scores.items(), key=lambda kv: kv[1], reverse=True)
sns.barplot(x=[str(r) for r, _ in cs_sorted], y=[s for _, s in cs_sorted],
            hue=[str(r) for r, _ in cs_sorted], palette="viridis", legend=False, ax=ax)
ax.set_xlabel("Residue")
ax.set_ylabel("Contact score")
ax.set_title("Contact scores by residue")
fig.tight_layout()
fig.savefig(PLOTS_DIR / "contact_scores.png", dpi=150)
plt.close(fig)

check(1, "len(ALL_COMBINATIONS) > 0 and <= 2048",
      0 < len(ALL_COMBINATIONS) <= 2048, actual=len(ALL_COMBINATIONS))
check(2, "FORCE_INCLUDE_RESIDUES present in every combination",
      all(all(r in c for r in FORCE_INCLUDE_RESIDUES) for c in ALL_COMBINATIONS))

# ---------------------------------------------------------------------------
# Cell 4 logic: mock 15 previous trial results
# ---------------------------------------------------------------------------
print("\n--- Cell 4: mock previous trial results ---")

rng = np.random.default_rng(42)
seen = set()
mock_trials = []
while len(mock_trials) < 15:
    combo = ALL_COMBINATIONS[rng.integers(len(ALL_COMBINATIONS))]
    if combo in seen:
        continue
    seen.add(combo)
    score = 15.0
    if 248 in combo:
        score -= 3.0
    if 312 in combo:
        score -= 2.0
    if 119 in combo:
        score -= 1.5
    score += rng.normal(0, 1.5)
    score = float(np.clip(score, 6.0, 20.0))
    mock_trials.append({
        "trial_id": len(mock_trials) + 1,
        "combination": list(combo),
        "combination_str": ",".join(map(str, combo)),
        "best_pae": score,
        "n_passed": int(rng.integers(0, 5)),
        "n_total": 5,
    })

trials_df = pd.DataFrame(mock_trials)
X_observed = np.array([encode(tuple(c)) for c in trials_df["combination"]])
y_observed = trials_df["best_pae"].to_numpy()

print(f"Trials completed: {len(trials_df)}")
print(f"Best pAE so far: {y_observed.min():.2f} (combination: "
      f"{trials_df.loc[trials_df['best_pae'].idxmin(), 'combination']})")
print(f"Worst pAE so far: {y_observed.max():.2f}")
print(f"Mean pAE: {y_observed.mean():.2f}")
print(trials_df.sort_values("best_pae")[["trial_id", "combination_str", "best_pae"]]
      .to_string(index=False))

check(3, "X_observed.shape == (15, 11)", X_observed.shape == (15, 11), actual=X_observed.shape)
check(4, "y_observed.shape == (15,)", y_observed.shape == (15,), actual=y_observed.shape)

# ---------------------------------------------------------------------------
# Cell 6 logic: fit Gaussian Process surrogate model
# ---------------------------------------------------------------------------
print("\n--- Cell 6: fit GP surrogate ---")

if len(X_observed) < 3:
    raise ValueError("Need at least 3 observations to fit a GP.")

kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5) + WhiteKernel(noise_level=0.1)
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5,
                               normalize_y=True, random_state=0)
gp.fit(X_observed, y_observed)

y_pred, y_std = gp.predict(ALL_COMBINATIONS_ENCODED, return_std=True)

tested_set = set(trials_df["combination_str"])
combo_strs = [",".join(map(str, c)) for c in ALL_COMBINATIONS]
already_tested = np.array([s in tested_set for s in combo_strs])

gp_df = pd.DataFrame({
    "combination": ALL_COMBINATIONS,
    "combination_str": combo_strs,
    "predicted_pae": y_pred,
    "uncertainty": y_std,
    "already_tested": already_tested,
    "contact_score_sum": COMBO_CONTACT_SCORES,
})

train_r2 = gp.score(X_observed, y_observed)
print("GP FIT SUMMARY")
print(f"Kernel: {gp.kernel_}")
print(f"Training points: {len(X_observed)}")
print(f"Training R^2: {train_r2:.2f}")
print("\nTop 10 predicted combinations (lowest predicted pAE):")
print(gp_df.sort_values("predicted_pae").head(10)
      [["combination_str", "predicted_pae", "uncertainty"]].to_string(index=False))

check(5, "GP fitted without error (gp.kernel_ is not None)", gp.kernel_ is not None)

fig, ax = plt.subplots(figsize=(7, 5))
untested_plot = gp_df[~gp_df["already_tested"]]
sns.scatterplot(data=untested_plot, x="predicted_pae", y="uncertainty",
                 hue="contact_score_sum", palette="viridis", ax=ax)
ax.set_title("GP landscape (untested combinations)")
fig.tight_layout()
fig.savefig(PLOTS_DIR / "gp_landscape.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------------
# Cell 7 logic: acquisition function and next batch selection
# ---------------------------------------------------------------------------
print("\n--- Cell 7: acquisition function & next batch ---")

untested_df = gp_df[~gp_df["already_tested"]].copy()
y_best = y_observed.min()
yp = untested_df["predicted_pae"].to_numpy()
ys = untested_df["uncertainty"].to_numpy()

if ACQUISITION_FUNC == "Expected Improvement":
    z = (y_best - yp - XI) / (ys + 1e-9)
    acq = (y_best - yp - XI) * norm.cdf(z) + ys * norm.pdf(z)
    acq = np.maximum(acq, 0)
elif ACQUISITION_FUNC == "Upper Confidence Bound":
    kappa = 2.0
    acq = -(yp - kappa * ys)
else:  # Thompson Sampling
    y_sample = rng.normal(yp, ys)
    acq = -y_sample

if USE_CONTACT_PRIOR:
    acq = acq * (0.5 + 0.5 * untested_df["contact_score_sum"].to_numpy())

untested_df["acquisition_score"] = acq
untested_df = untested_df.sort_values("acquisition_score", ascending=False)
NEXT_BATCH = list(untested_df.head(BATCH_SIZE)["combination"])

print(f"ACQUISITION FUNCTION: {ACQUISITION_FUNC} (xi={XI})")
print(f"Contact prior: {'enabled' if USE_CONTACT_PRIOR else 'disabled'}")
print("\nNEXT BATCH — top combinations selected:")
print(untested_df.head(BATCH_SIZE)
      [["combination_str", "predicted_pae", "uncertainty", "acquisition_score",
        "contact_score_sum"]].to_string(index=False))

check(6, "len(NEXT_BATCH) == BATCH_SIZE", len(NEXT_BATCH) == BATCH_SIZE, actual=len(NEXT_BATCH))

next_batch_strs = {",".join(map(str, c)) for c in NEXT_BATCH}
check(7, "No NEXT_BATCH combination already in mock trials",
      next_batch_strs.isdisjoint(tested_set))

contact_lookup = dict(zip(gp_df["combination_str"], gp_df["contact_score_sum"]))
check(10, "contact_score_sum > 0 for all NEXT_BATCH combinations",
      all(contact_lookup[s] > 0 for s in next_batch_strs))

fig, ax = plt.subplots(figsize=(9, 6))
top20 = untested_df.head(20)
sns.barplot(data=top20, x="acquisition_score", y="combination_str",
            hue="predicted_pae", dodge=False, palette="viridis", ax=ax, legend=False)
ax.set_title("Top 20 acquisition scores")
fig.tight_layout()
fig.savefig(PLOTS_DIR / "acquisition_scores.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------------
# Cell 8 logic: campaign progress dashboard (2x3 seaborn figure)
# ---------------------------------------------------------------------------
print("\n--- Cell 8: campaign progress dashboard ---")

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Plot 1: learning curve
running_min = trials_df["best_pae"].cummin()
sns.lineplot(x=trials_df["trial_id"], y=running_min, ax=axes[0, 0], label="Running best")
sns.scatterplot(x=trials_df["trial_id"], y=trials_df["best_pae"], ax=axes[0, 0], label="Trial score")
best_idx = trials_df["best_pae"].idxmin()
axes[0, 0].scatter(trials_df.loc[best_idx, "trial_id"], trials_df.loc[best_idx, "best_pae"],
                    marker="*", s=300, color="red", zorder=5, label="Current best")
axes[0, 0].set_title("Learning curve")
axes[0, 0].set_xlabel("Trial number")
axes[0, 0].set_ylabel("best_pae")
axes[0, 0].legend(fontsize=8)

# Plot 2: score distribution
sns.histplot(trials_df["best_pae"], kde=True, ax=axes[0, 1])
axes[0, 1].axvline(y_observed.min(), color="red", linestyle="--", label="Current best")
axes[0, 1].axvline(10.0, color="black", linestyle=":", label="Target threshold")
axes[0, 1].set_title("Score distribution")
axes[0, 1].legend(fontsize=8)

# Plot 3: residue importance (top 25% of trials)
top_quartile_cut = trials_df["best_pae"].quantile(0.25)
top_trials = trials_df[trials_df["best_pae"] <= top_quartile_cut]
importance = {r: 0 for r in candidate_pool}
for combo in top_trials["combination"]:
    for r in combo:
        if r in importance:
            importance[r] += 1
imp_df = pd.DataFrame({"residue": list(importance.keys()), "count": list(importance.values())})
imp_df["contact_score"] = imp_df["residue"].map(contact_scores)
imp_df = imp_df.sort_values("count", ascending=False)
sns.barplot(data=imp_df, x="residue", y="count", hue="contact_score",
            dodge=False, palette="viridis", ax=axes[0, 2], legend=False)
axes[0, 2].set_title("Residue importance (top 25% trials)")

# Plot 4: pairwise co-occurrence heatmap (top 25% of trials)
co_matrix = pd.DataFrame(0, index=candidate_pool, columns=candidate_pool)
for combo in top_trials["combination"]:
    for r1, r2 in itertools.combinations(combo, 2):
        co_matrix.loc[r1, r2] += 1
        co_matrix.loc[r2, r1] += 1
sns.heatmap(co_matrix, annot=True, fmt="d", cmap="viridis", ax=axes[1, 0], cbar=False)
axes[1, 0].set_title("Pairwise co-occurrence (top 25%)")

# Plot 5: score vs combination size
trials_df["n_residues"] = trials_df["combination"].apply(len)
sns.boxplot(data=trials_df, x="n_residues", y="best_pae", ax=axes[1, 1])
sns.stripplot(data=trials_df, x="n_residues", y="best_pae", color="black", alpha=0.5, ax=axes[1, 1])
axes[1, 1].set_title("Score vs combination size")

# Plot 6: GP predicted vs actual for tested combinations
gp_predicted_lookup = gp_df.set_index("combination_str")["predicted_pae"]
trials_df["gp_predicted"] = trials_df["combination_str"].map(gp_predicted_lookup)
sns.scatterplot(data=trials_df, x="gp_predicted", y="best_pae", ax=axes[1, 2])
lo = min(trials_df["gp_predicted"].min(), trials_df["best_pae"].min())
hi = max(trials_df["gp_predicted"].max(), trials_df["best_pae"].max())
axes[1, 2].plot([lo, hi], [lo, hi], color="red", linestyle="--")
axes[1, 2].set_title("GP predicted vs actual")

fig.tight_layout()
dashboard_path = PLOTS_DIR / "campaign_dashboard.png"
fig.savefig(dashboard_path, dpi=300)
plt.close(fig)

print(f"Dashboard saved: {dashboard_path}")
print(f"Current best: {y_observed.min():.2f} (trial {trials_df.loc[best_idx, 'trial_id']}, "
      f"combination {trials_df.loc[best_idx, 'combination']})")

n_panels = axes.shape[0] * axes.shape[1]
check(8, "All 6 dashboard plots generated and saved",
      dashboard_path.exists() and n_panels == 6, actual=(dashboard_path.exists(), n_panels))

# ---------------------------------------------------------------------------
# Cell 10 logic: check convergence
# ---------------------------------------------------------------------------
print("\n--- Cell 10: convergence check ---")

if len(y_observed) >= CONVERGENCE_PATIENCE:
    recent_best = min(y_observed[-CONVERGENCE_PATIENCE:])
    overall_best = min(y_observed)
    improvement = recent_best - overall_best
    CONVERGED = improvement < CONVERGENCE_THRESHOLD
else:
    CONVERGED = False

print("CONVERGENCE CHECK")
print(f"Trials completed: {len(trials_df)} / (budget not tracked in dry run)")
print(f"Best pAE: {y_observed.min():.2f}")
print(f"Converged: {'YES' if CONVERGED else 'NO'}")

check(9, "Convergence check ran (CONVERGED is defined)", isinstance(CONVERGED, (bool, np.bool_)))

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("FINAL REPORT")
print("=" * 70)
all_pass = all(v == "PASS" for v in results.values())
for n in sorted(results):
    print(f"  CHECK {n}: {results[n]}")
print(f"\nOverall: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
print(f"Plots written to: {PLOTS_DIR}")

if not all_pass:
    raise SystemExit(1)
"""Modeling pipeline (development set only). Every data-dependent step is fitted on the training part of each fold.

Schemes
  cv   : 5 repeats x StratifiedGroupKFold(5), y = dx|cohort, groups = split_group  (identical folds for every config)
  loso : leave-one-cohort-out over the task's cohorts (held-out cohort self-standardised, label-free)
Every (config, task, scheme, fold, seed) writes one JSON under results/runs/ and is skipped if it already exists.
"""
import json, os, time, hashlib, warnings
import numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.feature_selection import f_classif
from common import DATA, RES, write_json, versions
from evalstats import metrics

warnings.filterwarnings("ignore")
RUNS = RES / "runs"
AUDIT = json.load(open(DATA.parent / "DATA_AUDIT.json"))
TASKS = {k: {"cohorts": v["cohorts"], "classes": v["classes"]} for k, v in AUDIT["tasks"].items()}
TASK_ORDER = ["A_SCZ_vs_CTL", "A_BD_vs_CTL", "A_MDD_vs_CTL", "B_multiclass"]
CS = np.logspace(-5, 0, 6)
_VERS = None

def vers():
    global _VERS
    if _VERS is None:
        f = RES / "versions.json"
        if f.exists(): _VERS = json.load(open(f))
        else: _VERS = versions(); write_json(_VERS, f)
    return _VERS

# ------------------------------------------------------------------ data
_X = {}
def load_dev():
    if "dev" not in _X:
        X = pd.read_pickle(DATA / "features_dev.pkl")
        meta = pd.read_csv(DATA / "split_locked.tsv", sep="\t", index_col=0, keep_default_na=False).loc[X.index]
        meta["affected"] = meta.affected.astype(str) == "True"
        _X["dev"] = (X, meta)
    return _X["dev"]

def task_data(task, population="primary", drop_ids=(), exclude_cohorts=()):
    """population: 'primary' (affected removed) or 'all' (all eligible)."""
    X, meta = load_dev()
    t = TASKS[task]
    m = meta.cohort.isin(t["cohorts"]) & meta.dx.isin(t["classes"]) & ~meta.cohort.isin(exclude_cohorts)
    if population == "primary": m &= ~meta.affected
    m &= ~meta.index.isin(list(drop_ids))
    return X[m], meta[m]

def cv_folds(meta, n_rep=5, k=5):
    y = (meta.dx + "|" + meta.cohort).values
    out = []
    for r in range(n_rep):
        for i, (tr, te) in enumerate(StratifiedGroupKFold(k, shuffle=True, random_state=r).split(meta, y, meta.split_group)):
            out.append((f"r{r}f{i}", tr, te))
    return out

def loso_folds(meta):
    return [(f"loso_{c}", np.where(meta.cohort.values != c)[0], np.where(meta.cohort.values == c)[0])
            for c in sorted(meta.cohort.unique())]

# ------------------------------------------------------------------ batch adjustment (fit on train)
def batch_std(Xtr, btr, Xte, bte):
    Xtr, Xte = Xtr.copy(), Xte.copy()
    for b in np.unique(np.r_[btr, bte]):
        itr, ite = btr == b, bte == b
        if itr.sum() >= 2:
            mu, sd = Xtr[itr].mean(0), Xtr[itr].std(0) + 1e-6
            Xtr[itr] = (Xtr[itr] - mu) / sd
        else:                                  # unseen batch (LOSO): label-free self-standardisation
            mu, sd = Xte[ite].mean(0), Xte[ite].std(0) + 1e-6
        if ite.any(): Xte[ite] = (Xte[ite] - mu) / sd
    return Xtr, Xte

def combat(Xtr, btr, ytr, Xte, bte):
    """Add-on ComBat: parametric EB fitted on TRAINING donors with diagnosis as protected covariate (neuroCombat);
    held-out donors of seen batches adjusted with the training batch estimates and the training-average covariate
    term (no held-out labels). Unseen batches: label-free self-standardisation onto the training scale."""
    from neuroCombat import neuroCombat
    import contextlib, io
    cov = pd.DataFrame({"batch": btr, "dx": ytr})
    with contextlib.redirect_stdout(io.StringIO()):
        res = neuroCombat(dat=Xtr.T, covars=cov, batch_col="batch", categorical_cols=["dx"])
    est = res["estimates"]; Xtr_c = res["data"].T
    levels = [str(x) for x in est["batches"]]
    var_pooled = est["var.pooled"].ravel(); stand = est["stand.mean"][:, 0] + est["mod.mean"].mean(axis=1)
    Xte_c = Xte.copy()
    for b in np.unique(bte):
        ite = bte == b
        if str(b) in levels:
            j = levels.index(str(b))
            z = (Xte[ite] - stand) / np.sqrt(var_pooled)
            z = (z - est["gamma.star"][j]) / np.sqrt(est["delta.star"][j])
            Xte_c[ite] = z * np.sqrt(var_pooled) + stand
        else:
            mu, sd = Xte[ite].mean(0), Xte[ite].std(0) + 1e-6
            Xte_c[ite] = (Xte[ite] - mu) / sd * np.sqrt(var_pooled) + stand
    # common scale for the classifier (fitted on training)
    mu, sd = Xtr_c.mean(0), Xtr_c.std(0) + 1e-6
    return (Xtr_c - mu) / sd, (Xte_c - mu) / sd

# ------------------------------------------------------------------ graphs / gene sets
_G = {}
def string_adj():
    if "string" not in _G:
        e = np.load(DATA / "graph_string.npz"); _G["string"] = (e["edges"], int(e["n"]))
    return _G["string"]

def norm_adj(edges, n):
    A = sp.coo_matrix((np.ones(2 * len(edges)), (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])), shape=(n, n)).tocsr()
    A.data[:] = 1.0
    d = np.asarray(A.sum(1)).ravel(); d[d == 0] = 1
    Dm = sp.diags(1 / np.sqrt(d))
    return (Dm @ A @ Dm).tocsr()

def random_edges(n_edges, n, seed):
    """Erdos-Renyi graph with exactly n_edges distinct undirected edges on the same node set."""
    rng = np.random.default_rng(seed); got = set()
    while len(got) < n_edges:
        a = rng.integers(0, n, 2 * (n_edges - len(got))); b = rng.integers(0, n, len(a))
        for x, y in zip(a, b):
            if x != y:
                got.add((min(x, y), max(x, y)))
                if len(got) == n_edges: break
    return np.array(sorted(got))

def coexpr_edges(Xtr, k=10):
    """kNN co-expression graph from TRAINING donors: each gene linked to its k genes of largest |Pearson r|; symmetrised."""
    Z = (Xtr - Xtr.mean(0)) / (Xtr.std(0) + 1e-9)
    n = Z.shape[1]; nb = np.zeros((n, k), dtype=np.int64)
    for s in range(0, n, 2000):
        C = np.abs(Z[:, s:s + 2000].T @ Z) / Z.shape[0]
        for i in range(C.shape[0]): C[i, s + i] = -1
        nb[s:s + 2000] = np.argpartition(-C, k, axis=1)[:, :k]
    e = np.c_[np.repeat(np.arange(n), k), nb.ravel()]
    e = np.unique(np.sort(e, axis=1), axis=0)
    return e[e[:, 0] != e[:, 1]]

def smooth(Xtr, Xte, An, alpha=0.5):
    Ytr = (1 - alpha) * Xtr + alpha * np.asarray((An @ Xtr.T).T)
    Yte = (1 - alpha) * Xte + alpha * np.asarray((An @ Xte.T).T)
    mu, sd = Ytr.mean(0), Ytr.std(0) + 1e-6          # re-standardise on training (no-op for no-graph)
    return (Ytr - mu) / sd, (Yte - mu) / sd

_PW = {}
def reactome_sets(genes):
    if "sets" not in _PW:
        r = pd.read_csv(DATA / "knowledge" / "NCBI2Reactome.txt", sep="\t", header=None, dtype=str,
                        names=["entrez", "pid", "url", "name", "evidence", "species"])
        r = r[r.species == "Homo sapiens"]
        gi = {g: i for i, g in enumerate(genes)}
        sets = {}
        for pid, grp in r.groupby("pid"):
            idx = sorted({gi[g] for g in grp.entrez if g in gi})
            if 10 <= len(idx) <= 300: sets[pid] = np.array(idx)
        _PW["sets"] = sets
    return _PW["sets"]

def celltype_sets(genes):
    """BrainInABlender (Hagenauer 2018) marker database: genes per CellType_Primary, human symbols mapped to Entrez;
    genes listed for >1 primary cell type dropped; >= 5 genes per set. Simplification vs Sir_UnMixALot: equal weight per
    gene (not per-publication averaging)."""
    if "ct" not in _PW:
        from common import symbols
        mk = pd.read_csv(DATA / "knowledge" / "BrainInABlender_markers.tsv", sep="\t", dtype=str)
        s2e = symbols().reset_index().dropna(subset=["SYMBOL"]).drop_duplicates("SYMBOL").set_index("SYMBOL")["ENTREZID"]
        mk["entrez"] = mk.GeneSymbol_Human.map(s2e)
        mk = mk.dropna(subset=["entrez"]).drop_duplicates(["entrez", "CellType_Primary"])
        multi = mk.groupby("entrez").CellType_Primary.nunique(); mk = mk[mk.entrez.map(multi) == 1]
        gi = {g: i for i, g in enumerate(genes)}
        sets = {ct: np.array(sorted({gi[g] for g in grp.entrez if g in gi})) for ct, grp in mk.groupby("CellType_Primary")}
        _PW["ct"] = {k: v for k, v in sets.items() if len(v) >= 5}
    return _PW["ct"]

def set_scores(Xtr, Xte, sets):
    M = np.stack([Xtr[:, s].mean(1) for s in sets.values()], 1); N = np.stack([Xte[:, s].mean(1) for s in sets.values()], 1)
    mu, sd = M.mean(0), M.std(0) + 1e-6
    return (M - mu) / sd, (N - mu) / sd

def random_sets(sets, n_genes, seed):
    rng = np.random.default_rng(seed)
    return {k: rng.choice(n_genes, len(v), replace=False) for k, v in sets.items()}

# ------------------------------------------------------------------ models
def fit_predict(model, Xtr, ytr, Xte, gtr, seed):
    classes = np.unique(ytr)
    if model == "lr":
        # Exact speed-up: the L2-penalised solution lies in the span of the training rows (representer theorem), so
        # fitting on R = X V (V = right singular vectors of the outer-training X) gives identical predictions, also for
        # the inner-CV fits (their rows are a subset). Verified: max |dP| = 3e-16 at the selected C (LOG.md).
        if Xtr.shape[1] > Xtr.shape[0]:
            _, _, Vt = np.linalg.svd(Xtr, full_matrices=False); Xtr, Xte = Xtr @ Vt.T, Xte @ Vt.T
        inner = list(StratifiedGroupKFold(3, shuffle=True, random_state=seed).split(Xtr, ytr, gtr))
        clf = LogisticRegressionCV(Cs=CS, cv=inner, class_weight="balanced", scoring="neg_log_loss", max_iter=10000, tol=1e-8,
                                   l1_ratios=(0,), use_legacy_attributes=False)
        clf.fit(Xtr, ytr); info = {"C": float(np.atleast_1d(clf.C_)[0])}
    elif model == "rf":
        clf = RandomForestClassifier(500, max_features="sqrt", class_weight="balanced", n_jobs=1, random_state=seed)
        clf.fit(Xtr, ytr); info = {}
    else:
        raise ValueError(model)
    P = clf.predict_proba(Xte)
    return clf.classes_, P, info

# ------------------------------------------------------------------ one fold
def run_fold(config, spec, task, scheme, fold_id, tr, te, X, meta, seed, tag=""):
    out_f = RUNS / config / task / f"{scheme}_{fold_id}_s{seed}{tag}.json"
    if out_f.exists(): return json.load(open(out_f))
    t0 = time.time()
    Xtr, Xte = X.values[tr], X.values[te]
    mtr, mte = meta.iloc[tr], meta.iloc[te]
    ytr, yte = mtr.dx.values, mte.dx.values
    info = {}
    if spec.get("model") == "study_only":
        cohorts = sorted(meta.cohort.unique())
        oh = lambda m: np.stack([(m.cohort.values == c).astype(float) for c in cohorts], 1)
        clf = LogisticRegression(C=1e4, class_weight="balanced", max_iter=5000).fit(oh(mtr), ytr)
        classes, P = clf.classes_, clf.predict_proba(oh(mte))
    else:
        # 1. batch adjustment
        if spec.get("adjust", "bstd") == "bstd": Xtr, Xte = batch_std(Xtr, mtr.batch.values, Xte, mte.batch.values)
        elif spec["adjust"] == "combat_infold": Xtr, Xte = combat(Xtr, mtr.batch.values, ytr, Xte, mte.batch.values)
        elif spec["adjust"] == "combat_alldata":   # LEAKAGE DEMO ONLY: ComBat fitted on train+test with ALL labels
            Xall = np.r_[Xtr, Xte]; b = np.r_[mtr.batch.values, mte.batch.values]; y = np.r_[ytr, yte]
            from neuroCombat import neuroCombat
            import contextlib, io
            with contextlib.redirect_stdout(io.StringIO()):
                Z = neuroCombat(dat=Xall.T, covars=pd.DataFrame({"batch": b, "dx": y}), batch_col="batch", categorical_cols=["dx"])["data"].T
            mu, sd = Z[:len(tr)].mean(0), Z[:len(tr)].std(0) + 1e-6
            Xtr, Xte = (Z[:len(tr)] - mu) / sd, (Z[len(tr):] - mu) / sd
        elif spec["adjust"] == "none":
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6; Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
        # 2. feature transform
        tf = spec.get("transform", "none")
        if tf in ("string", "string_random", "coexpr", "coexpr_random"):
            n = Xtr.shape[1]
            if tf.startswith("string"): E, _ = string_adj()
            else: E = coexpr_edges(Xtr, k=spec.get("k", 10))
            if tf.endswith("random"): E = random_edges(len(E), n, seed=hash_seed(config, task, fold_id, seed))
            info["n_edges"] = int(len(E))
            Xtr, Xte = smooth(Xtr, Xte, norm_adj(E, n), spec.get("alpha", 0.5))
        elif tf in ("reactome", "reactome_random"):
            sets = reactome_sets(list(X.columns))
            if tf.endswith("random"): sets = random_sets(sets, X.shape[1], hash_seed(config, task, fold_id, seed))
            info["n_sets"] = len(sets)
            Xtr, Xte = set_scores(Xtr, Xte, sets)
        elif tf in ("celltype", "celltype_random"):
            sets = celltype_sets(list(X.columns))
            if tf.endswith("random"): sets = random_sets(sets, X.shape[1], hash_seed(config, task, fold_id, seed))
            info["n_sets"] = len(sets)
            Xtr, Xte = set_scores(Xtr, Xte, sets)
        if spec.get("filter_k"):
            F, _ = f_classif(Xtr, ytr); keep = np.argsort(-np.nan_to_num(F))[:spec["filter_k"]]
            Xtr, Xte = Xtr[:, keep], Xte[:, keep]
        classes, P, inf = fit_predict(spec.get("model", "lr"), Xtr, ytr, Xte, mtr.split_group.values, seed)
        info.update(inf)
    pred = classes[P.argmax(1)]
    res = {"config": config, "task": task, "scheme": scheme, "fold": fold_id, "seed": seed, "tag": tag, "spec": spec,
           "n_train": int(len(tr)), "n_test": int(len(te)),
           "train_class_counts": {k: int(v) for k, v in pd.Series(ytr).value_counts().items()},
           "test_class_counts": {k: int(v) for k, v in pd.Series(yte).value_counts().items()},
           "test_cohorts": sorted(set(mte.cohort)), "classes": [str(c) for c in classes],
           "donors": list(mte.index), "y": list(yte), "pred": list(pred), "proba": np.round(P, 5).tolist(),
           "metrics": metrics(yte, pred, P, classes), "info": info, "seconds": round(time.time() - t0, 2), "versions": vers()}
    write_json(res, out_f)
    return res

def hash_seed(*a):
    return int(hashlib.md5("|".join(map(str, a)).encode()).hexdigest()[:8], 16)

def run_config(config, spec, tasks=TASK_ORDER, schemes=("cv", "loso"), seeds=(0,), population="primary", n_jobs=4, tag="", **kw):
    from joblib import Parallel, delayed
    jobs = []
    for task in tasks:
        X, meta = task_data(task, population=population, **kw)
        for scheme in schemes:
            folds = cv_folds(meta) if scheme == "cv" else loso_folds(meta)
            for fid, tr, te in folds:
                for s in seeds:
                    jobs.append((task, scheme, fid, tr, te, X, meta, s))
    Parallel(n_jobs=n_jobs, backend="loky")(delayed(run_fold)(config, spec, t, sc, f, tr, te, X, meta, s, tag) for t, sc, f, tr, te, X, meta, s in jobs)
    return len(jobs)

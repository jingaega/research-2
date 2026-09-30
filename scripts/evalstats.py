"""Evaluation statistics: chance levels, Nadeau-Bengio corrected resampled t-test, bootstrap CIs."""
import numpy as np
from scipy import stats
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix, precision_recall_fscore_support

def chance_levels(y_eval, y_train=None, n_sim=5000, seed=0):
    """Expected macro-F1 on the evaluation labels of: (i) uniform random guessing, (ii) random guessing
    with training-set class priors, (iii) always predicting the training majority class. Simulated exactly
    on the actual evaluation label vector."""
    y_eval = np.asarray(y_eval); classes = np.unique(y_eval if y_train is None else np.concatenate([y_eval, y_train]))
    pri = np.array([(np.asarray(y_train if y_train is not None else y_eval) == c).mean() for c in classes])
    rng = np.random.default_rng(seed)
    uni = [f1_score(y_eval, rng.choice(classes, len(y_eval)), average="macro", labels=classes, zero_division=0) for _ in range(n_sim)]
    pr = [f1_score(y_eval, rng.choice(classes, len(y_eval), p=pri), average="macro", labels=classes, zero_division=0) for _ in range(n_sim)]
    maj = f1_score(y_eval, np.full(len(y_eval), classes[pri.argmax()]), average="macro", labels=classes, zero_division=0)
    return {"uniform_random_macroF1": float(np.mean(uni)), "prior_random_macroF1": float(np.mean(pr)),
            "majority_class_macroF1": float(maj), "class_priors": dict(zip(map(str, classes), pri.round(4).tolist())),
            "mcc_chance": 0.0, "auc_chance": 0.5}

def metrics(y, pred, proba=None, classes=None):
    classes = np.unique(y) if classes is None else np.asarray(classes)
    p, r, f, sup = precision_recall_fscore_support(y, pred, labels=classes, zero_division=0)
    out = {"macro_f1": float(f1_score(y, pred, average="macro", labels=classes, zero_division=0)),
           "mcc": float(matthews_corrcoef(y, pred)),
           "per_class": {str(c): {"precision": float(p[i]), "recall": float(r[i]), "f1": float(f[i]), "n": int(sup[i])} for i, c in enumerate(classes)},
           "confusion_matrix": {"labels": [str(c) for c in classes], "matrix": confusion_matrix(y, pred, labels=classes).tolist()}}
    if proba is not None:
        try:
            if len(classes) == 2: out["auc"] = float(roc_auc_score((np.asarray(y) == classes[1]).astype(int), proba[:, 1]))
            else: out["auc"] = float(roc_auc_score(y, proba, multi_class="ovr", average="macro", labels=classes))
        except ValueError: out["auc"] = None
    return out

def nadeau_bengio(diffs, n_train, n_test):
    """Corrected resampled t-test (Nadeau & Bengio 2003; Bouckaert & Frank 2004 form for r x k CV) on the
    k*r fold-level differences. Assumes the correlation between fold-level estimates induced by overlapping
    training sets is approximated by n_test/n_train. Returns corrected and uncorrected two-sided p."""
    d = np.asarray(diffs, float); J = len(d); mu = d.mean(); var = d.var(ddof=1)
    if var == 0: return {"mean_diff": mu, "t": np.nan, "p_corrected": np.nan, "p_uncorrected": np.nan, "J": J}
    t_c = mu / np.sqrt((1 / J + n_test / n_train) * var); t_u = mu / np.sqrt(var / J)
    return {"mean_diff": float(mu), "sd_diff": float(np.sqrt(var)), "J": J, "df": J - 1, "t_corrected": float(t_c),
            "p_corrected": float(2 * stats.t.sf(abs(t_c), J - 1)), "t_uncorrected": float(t_u),
            "p_uncorrected": float(2 * stats.t.sf(abs(t_u), J - 1)), "n_train": float(n_train), "n_test": float(n_test),
            "assumption": "between-fold correlation approximated by n_test/n_train (Nadeau-Bengio heuristic)"}

def bootstrap_ci(y, pred, proba=None, n_boot=2000, seed=0, classes=None):
    """Percentile bootstrap (resampling donors of the locked test set) for macro-F1, MCC, AUC."""
    y = np.asarray(y); pred = np.asarray(pred); rng = np.random.default_rng(seed); n = len(y)
    classes = np.unique(y) if classes is None else classes
    vals = {"macro_f1": [], "mcc": [], "auc": []}
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if len(np.unique(y[i])) < len(classes): continue
        vals["macro_f1"].append(f1_score(y[i], pred[i], average="macro", labels=classes, zero_division=0))
        vals["mcc"].append(matthews_corrcoef(y[i], pred[i]))
        if proba is not None:
            try:
                vals["auc"].append(roc_auc_score((y[i] == classes[1]).astype(int), proba[i, 1]) if len(classes) == 2
                                   else roc_auc_score(y[i], proba[i], multi_class="ovr", average="macro", labels=classes))
            except ValueError: pass
    return {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v else None for k, v in vals.items()}

"""Rule 12 runtime measurement. Uses DEV donors only and RANDOMLY PERMUTED labels, and computes NO performance
metric: only wall-clock time is recorded (PLAN.md not yet approved, so nothing may be evaluated).
Output: results/timing.json"""
import time, json
import numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedGroupKFold
from common import DATA, RES, write_json, versions

X = pd.read_pickle(DATA / "features_dev.pkl")
meta = pd.read_csv(DATA / "split_locked.tsv", sep="\t", index_col=0, keep_default_na=False).loc[X.index]
m = meta.cohort.isin(["Stanley", "Pittsburgh", "Pritzker", "CharingCross", "Victoria"]) & meta.dx.isin(["SCZ", "CTL"])
X, meta = X[m], meta[m]
rng = np.random.default_rng(0); y = rng.permutation(meta.dx.values)   # permuted: no real signal can be seen
E = np.load(DATA / "graph_string.npz")["edges"]; n = X.shape[1]
A = sp.coo_matrix((np.ones(len(E) * 2), (np.r_[E[:, 0], E[:, 1]], np.r_[E[:, 1], E[:, 0]])), shape=(n, n)).tocsr()
d = np.asarray(A.sum(1)).ravel(); d[d == 0] = 1; Dm = sp.diags(1 / np.sqrt(d)); An = Dm @ A @ Dm

def batch_std(Xtr, btr, Xte, bte):
    Xtr, Xte = Xtr.copy(), Xte.copy()
    for b in np.unique(btr):
        mu = Xtr[btr == b].mean(0); sd = Xtr[btr == b].std(0) + 1e-6
        Xtr[btr == b] = (Xtr[btr == b] - mu) / sd
        if (bte == b).any(): Xte[bte == b] = (Xte[bte == b] - mu) / sd
    return Xtr, Xte

T = {}
for name, smooth in [("LR_no_graph", False), ("LR_graph_smoothed", True)]:
    t0 = time.time(); nf = 0
    for rep in range(5):
        for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=rep).split(X, meta.dx + meta.cohort, meta.split_group):
            Xtr, Xte = batch_std(X.values[tr], meta.batch.values[tr], X.values[te], meta.batch.values[te])
            if smooth: Xtr = 0.5 * Xtr + 0.5 * (An @ Xtr.T).T; Xte = 0.5 * Xte + 0.5 * (An @ Xte.T).T
            LogisticRegressionCV(Cs=np.logspace(-5, 0, 6), cv=3, class_weight="balanced", max_iter=2000, n_jobs=4).fit(Xtr, y[tr])
            nf += 1
    T[name] = {"seconds_25_folds": round(time.time() - t0, 1), "folds": nf}
    print(name, T[name])
# one GCN fold (2-layer graph conv over genes, hidden 8, mean-pool, 150 epochs) to extrapolate
import torch
torch.set_num_threads(4)
tr, te = next(StratifiedGroupKFold(5, shuffle=True, random_state=0).split(X, meta.dx + meta.cohort, meta.split_group))
Xtr, _ = batch_std(X.values[tr], meta.batch.values[tr], X.values[te], meta.batch.values[te])
Ahat = An + sp.eye(n); Ahat = Ahat.tocoo()
At = torch.sparse_coo_tensor(np.vstack([Ahat.row, Ahat.col]), Ahat.data.astype(np.float32), (n, n)).coalesce()
xt = torch.tensor(Xtr, dtype=torch.float32); yt = torch.tensor((y[tr] == "SCZ").astype(np.int64))
class GCN(torch.nn.Module):
    def __init__(s, h=8):
        super().__init__(); s.w1 = torch.nn.Linear(1, h); s.w2 = torch.nn.Linear(h, h); s.out = torch.nn.Linear(h, 2)
    def forward(s, x):                       # x: batch x genes
        B = x.shape[0]; H = x.T.unsqueeze(-1)                      # genes x batch x 1
        H = torch.relu(s.w1(torch.sparse.mm(At, H.reshape(n, -1)).reshape(n, B, 1)))
        H = torch.relu(s.w2(torch.sparse.mm(At, H.reshape(n, -1)).reshape(n, B, -1)))
        return s.out(H.mean(0))
net = GCN(); opt = torch.optim.Adam(net.parameters(), 1e-2, weight_decay=1e-4)
t0 = time.time()
for ep in range(150):
    for i in range(0, len(xt), 64):
        opt.zero_grad(); loss = torch.nn.functional.cross_entropy(net(xt[i:i + 64]), yt[i:i + 64]); loss.backward(); opt.step()
T["GCN_one_fold_150_epochs"] = {"seconds": round(time.time() - t0, 1)}
print(T["GCN_one_fold_150_epochs"])
write_json({"versions": versions(), "note": "dev donors only, permuted labels, no metrics computed", "task": "SCZ vs CTL (largest task)",
            "n_donors": int(len(X)), "n_genes": int(n), "timings": T}, RES / "timing.json")

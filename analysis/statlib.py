"""Minimal statistics helpers (numpy/scipy only) for the revision analyses."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats, optimize


def z(x):
    x = pd.to_numeric(x, errors="coerce")
    return (x - x.mean()) / x.std(ddof=1)


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1)
    out[ok] = adj
    return out


def design(df, terms, cat=()):
    X = pd.DataFrame(index=df.index)
    X["const"] = 1.0
    for t in terms:
        X[t] = pd.to_numeric(df[t], errors="coerce")
    for c in cat:
        dm = pd.get_dummies(df[c].astype("category"), prefix=c, drop_first=True, dtype=float)
        dm[df[c].isna()] = np.nan
        X = pd.concat([X, dm], axis=1)
    return X


def ols(y, X):
    """y: Series, X: DataFrame incl const. complete cases. t-based inference."""
    d = pd.concat([y.rename("_y"), X], axis=1).dropna()
    yy = d["_y"].to_numpy(float)
    XX = d.drop(columns="_y")
    XX = XX.loc[:, XX.std(ddof=0).gt(0) | (XX.columns == "const")]
    Xm = XX.to_numpy(float)
    n, p = Xm.shape
    beta, *_ = np.linalg.lstsq(Xm, yy, rcond=None)
    resid = yy - Xm @ beta
    rss = float(resid @ resid)
    tss = float(((yy - yy.mean()) ** 2).sum())
    dfr = n - p
    s2 = rss / dfr
    cov = np.linalg.pinv(Xm.T @ Xm) * s2
    se = np.sqrt(np.diag(cov))
    t = beta / se
    pval = 2 * stats.t.sf(np.abs(t), dfr)
    tcrit = stats.t.ppf(0.975, dfr)
    r2 = 1 - rss / tss
    adj = 1 - (1 - r2) * (n - 1) / dfr
    ll = -n / 2 * (np.log(2 * np.pi * rss / n) + 1)
    res = pd.DataFrame({"beta": beta, "se": se, "lo": beta - tcrit * se, "hi": beta + tcrit * se, "p": pval}, index=XX.columns)
    return dict(coef=res, n=n, p=p, rss=rss, r2=r2, adj_r2=adj, aic=-2 * ll + 2 * (p + 1), bic=-2 * ll + np.log(n) * (p + 1),
                cov=pd.DataFrame(cov, index=XX.columns, columns=XX.columns), dfr=dfr, idx=d.index)


def lincom(fit, w: dict):
    """Linear combination of coefficients with t CI."""
    names = list(w)
    v = np.array([w[k] for k in names])
    b = fit["coef"].loc[names, "beta"].to_numpy()
    C = fit["cov"].loc[names, names].to_numpy()
    est = float(v @ b)
    se = float(np.sqrt(v @ C @ v))
    tcrit = stats.t.ppf(0.975, fit["dfr"])
    p = 2 * stats.t.sf(abs(est / se), fit["dfr"])
    return est, est - tcrit * se, est + tcrit * se, p


def logit(y, X, maxiter=100):
    d = pd.concat([y.rename("_y"), X], axis=1).dropna()
    yy = d["_y"].to_numpy(float)
    Xm = d.drop(columns="_y").to_numpy(float)
    b = np.zeros(Xm.shape[1])
    for _ in range(maxiter):
        eta = Xm @ b
        mu = 1 / (1 + np.exp(-eta))
        W = mu * (1 - mu)
        g = Xm.T @ (yy - mu)
        H = Xm.T @ (Xm * W[:, None])
        step = np.linalg.solve(H + 1e-9 * np.eye(len(b)), g)
        b = b + step
        if np.max(np.abs(step)) < 1e-9:
            break
    mu = 1 / (1 + np.exp(-(Xm @ b)))
    H = Xm.T @ (Xm * (mu * (1 - mu))[:, None])
    cov = np.linalg.inv(H)
    se = np.sqrt(np.diag(cov))
    ll = float(np.sum(yy * np.log(mu + 1e-15) + (1 - yy) * np.log(1 - mu + 1e-15)))
    cols = d.drop(columns="_y").columns
    res = pd.DataFrame({"beta": b, "se": se, "p": 2 * stats.norm.sf(np.abs(b / se))}, index=cols)
    return dict(coef=res, prob=mu, y=yy, ll=ll, n=len(yy), k=len(b), idx=d.index, X=Xm, b=b)


def auc(y, s):
    y = np.asarray(y); s = np.asarray(s)
    pos = s[y == 1]; neg = s[y == 0]
    r = stats.rankdata(np.concatenate([pos, neg]))
    return (r[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def delong(y, s1, s2):
    """DeLong test for two correlated AUCs."""
    y = np.asarray(y); s1 = np.asarray(s1); s2 = np.asarray(s2)
    pos = y == 1; neg = ~pos
    m, n = pos.sum(), neg.sum()

    def comps(s):
        X, Y = s[pos], s[neg]
        V10 = np.array([(np.sum(x > Y) + 0.5 * np.sum(x == Y)) / n for x in X])
        V01 = np.array([(np.sum(X > yv) + 0.5 * np.sum(X == yv)) / m for yv in Y])
        return V10.mean(), V10, V01

    a1, v10a, v01a = comps(s1)
    a2, v10b, v01b = comps(s2)
    S10 = np.cov(np.vstack([v10a, v10b]))
    S01 = np.cov(np.vstack([v01a, v01b]))
    S = S10 / m + S01 / n
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    zz = (a2 - a1) / np.sqrt(var)
    se1 = np.sqrt(S[0, 0]); se2 = np.sqrt(S[1, 1])
    return dict(auc1=a1, auc1_lo=a1 - 1.96 * se1, auc1_hi=a1 + 1.96 * se1,
                auc2=a2, auc2_lo=a2 - 1.96 * se2, auc2_hi=a2 + 1.96 * se2,
                diff=a2 - a1, diff_lo=a2 - a1 - 1.96 * np.sqrt(var), diff_hi=a2 - a1 + 1.96 * np.sqrt(var),
                p=2 * stats.norm.sf(abs(zz)))


def nri_idi(y, p1, p2):
    y = np.asarray(y); p1 = np.asarray(p1); p2 = np.asarray(p2)
    ev, ne = y == 1, y == 0
    up, dn = p2 > p1, p2 < p1
    nri_e = up[ev].mean() - dn[ev].mean()
    nri_n = dn[ne].mean() - up[ne].mean()
    idi = (p2[ev].mean() - p2[ne].mean()) - (p1[ev].mean() - p1[ne].mean())
    return nri_e + nri_n, nri_e, nri_n, idi


# ---------- Linear mixed model: random intercept + random slope (unstructured), REML ----------

def lmm_fit(y, X, Z, groups, reml=True):
    """y (n,), X (n,p) fixed design, Z (n,q) random design, groups (n,) labels.
    Returns fixed effects, cov (model-based), variance components."""
    y = np.asarray(y, float); X = np.asarray(X, float); Z = np.asarray(Z, float)
    g = pd.Series(groups).astype("category").cat.codes.to_numpy()
    G = g.max() + 1
    idx = [np.where(g == k)[0] for k in range(G)]
    q = Z.shape[1]
    n, p = X.shape

    def unpack(th):
        L = np.zeros((q, q))
        L[np.tril_indices(q)] = th
        return L

    def negll(th, ret=False):
        L = unpack(th)  # relative cholesky factor (D / sigma2 = L L')
        D = L @ L.T
        XtVX = np.zeros((p, p)); XtVy = np.zeros(p); ldet = 0.0
        blocks = []
        for ii in idx:
            Zi = Z[ii]; Xi = X[ii]; yi = y[ii]
            Vi = Zi @ D @ Zi.T + np.eye(len(ii))
            ci = np.linalg.cholesky(Vi)
            ldet += 2 * np.log(np.diag(ci)).sum()
            Vinv = np.linalg.inv(Vi)
            XtVX += Xi.T @ Vinv @ Xi
            XtVy += Xi.T @ Vinv @ yi
            blocks.append(Vinv)
        beta = np.linalg.solve(XtVX, XtVy)
        rss = 0.0
        for ii, Vinv in zip(idx, blocks):
            r = y[ii] - X[ii] @ beta
            rss += r @ Vinv @ r
        if reml:
            s2 = rss / (n - p)
            val = 0.5 * (ldet + (n - p) * np.log(s2) + np.linalg.slogdet(XtVX)[1] + (n - p))
        else:
            s2 = rss / n
            val = 0.5 * (ldet + n * np.log(s2) + n)
        if ret:
            return beta, s2, XtVX, D
        return val

    th0 = np.eye(q)[np.tril_indices(q)] * 0.5
    opt = optimize.minimize(negll, th0, method="L-BFGS-B")
    opt = optimize.minimize(negll, opt.x, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-7, fatol=1e-9))
    beta, s2, XtVX, D = negll(opt.x, ret=True)
    cov = np.linalg.inv(XtVX) * s2
    return dict(beta=beta, se=np.sqrt(np.diag(cov)), cov=cov, sigma2=s2, D=D * s2, n=n, groups=G, converged=opt.success,
                negll=opt.fun)

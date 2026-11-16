"""
EVITAR-TB protocol: planning numbers with the 30-day window period built in.
Primary-outcome events counted in (30 days, 24 months]; censoring (33%/yr)
runs from randomization; constant hazards; cumulative observed-event normal
approximation; B+MMS vs placebo at two-sided alpha 0.05 (fixed sequence,
B alone tested at 0.05 only after H1 is rejected). 3 arms, 1:1:1, N = 8,650.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

N, T, T0 = 8650, 2.0, 30 / 365.25

def p_obs(hr, risk=0.03, attr=0.33, t0=T0):
    lam, d = -math.log(1 - risk), -math.log(1 - attr)
    h = hr * lam
    return h / (h + d) * (math.exp(-(h + d) * t0) - math.exp(-(h + d) * T))

def power(hr, alpha=0.05, n=N, **k):
    p0, p1 = p_obs(1, **k), p_obs(hr, **k)
    se = math.sqrt((p0*(1-p0) + p1*(1-p1)) / (n / 3))
    return norm.cdf((p0 - p1) / se - norm.ppf(1 - alpha / 2))

def n_needed(hr, pw, alpha=0.05, raw=False, **k):
    p0, p1 = p_obs(1, **k), p_obs(hr, **k)
    z = norm.ppf(1 - alpha / 2) + norm.ppf(pw)
    x = 3 * z**2 * (p0*(1-p0) + p1*(1-p1)) / (p0 - p1)**2
    return x if raw else int(math.ceil(x / 50) * 50)

def ev(hr, **k):
    return N / 3 * p_obs(hr, **k)

def sim(hr_c, hr_b, reps=2_000_000, seed=108350):
    rng = np.random.default_rng(seed)
    n = N // 3; n0 = N - 2 * n
    p0, p1, p2 = p_obs(1), p_obs(hr_c), p_obs(hr_b)
    x0 = rng.binomial(n0, p0, reps) / n0
    x1 = rng.binomial(n, p1, reps) / n
    x2 = rng.binomial(n, p2, reps) / n
    def z(xa):
        se = np.sqrt(x0*(1-x0)/n0 + xa*(1-xa)/n)
        return np.divide(x0 - xa, se, out=np.zeros_like(se), where=se > 0)
    zc, zb, c = z(x1), z(x2), norm.ppf(0.975)
    r1 = (np.abs(zc) > c) & (zc > 0)
    r2 = r1 & (np.abs(zb) > c) & (zb > 0)
    return r1.mean(), r2.mean()

if __name__ == "__main__":
    hr80 = brentq(lambda h: power(h) - 0.8, 0.4, 0.95)
    print(f"detectable HR (80%): {hr80:.4f} = {1-hr80:.1%}; unrounded N at HR 0.656: {n_needed(0.656,0.8,raw=True):.0f}; "
          f"power at 0.656: {power(0.656):.4f}; at 0.65: {power(0.65):.4f}")
    print(f"events after day 30: placebo {ev(1):.1f}, HR0.656 {ev(0.656):.1f}, HR0.65 {ev(0.65):.1f}")
    print(f"total events both active HR0.65: {ev(1)+2*ev(0.65):.0f}; only combo effective: {2*ev(1)+ev(0.65):.0f}")
    old = brentq(lambda h: power(h, alpha=0.035) - 0.8, 0.4, 0.95)
    print(f"old design (B alone first, alpha 0.035) detectable with window: {1-old:.1%}")
    print("Table 20:")
    for red in (0.25, 0.30, 0.35, 0.40, 0.50):
        hr = round(1 - red, 2)
        print(f"  {red:.0%} {hr}: N80 {n_needed(hr,0.8):,}, N90 {n_needed(hr,0.9):,}, power {power(hr):.1%}")
    hr = 0.656
    print(f"  34.4% {hr}: N80 {n_needed(hr,0.8):,}, N90 {n_needed(hr,0.9):,}, power {power(hr):.1%}")
    print("Table 21 (HR 0.656):")
    for r in (0.02, 0.025, 0.03, 0.035, 0.04):
        print(f"  risk {r:.1%}: {power(hr, risk=r):.1%}, {ev(1, risk=r):.0f} / {ev(hr, risk=r):.0f}")
    for l in (0.15, 0.25, 0.33, 0.40, 0.50):
        print(f"  loss {l:.0%}: {power(hr, attr=l):.1%}, {ev(1, attr=l):.0f} / {ev(hr, attr=l):.0f}")
    print("Table 19 / synopsis (B alone power under fixed sequence):")
    for hc, hb in ((0.65, 0.65), (0.65, 0.70)):
        a, b = sim(hc, hb)
        print(f"  B+MMS {hc}, B {hb}: primary {a:.1%}, B alone {b:.1%}")

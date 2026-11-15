"""MODIS + VIIRS harmonization pipeline -> unified fire calendar + dashboard.
Usage:  python run.py              (reads FIRMS files from ./data)
        python run.py --demo       (synthetic data, only to test that everything runs)
"""
import argparse, glob, json, os, sys, webbrowser, zipfile
import numpy as np, pandas as pd

VCONF = {"l": .2, "low": .2, "n": .6, "nominal": .6, "h": .9, "high": .9}
CODE = {"MODIS_only": 1, "VIIRS_only": 2, "both": 3}

def read_files(folder):
    found = []
    for p in sorted(glob.glob(os.path.join(folder, "**", "*"), recursive=True)):
        low = p.lower()
        if low.endswith(".zip"):
            with zipfile.ZipFile(p) as z:
                for n in z.namelist():
                    if n.lower().endswith(".csv"):
                        found.append((os.path.basename(p), pd.read_csv(z.open(n))))
        elif low.endswith(".csv"):
            found.append((os.path.basename(p), pd.read_csv(p)))
    return [(n, d) for n, d in found if {"acq_date", "acq_time", "latitude"} <= set(d.columns)]

def normalize(d, sensor):
    t = d["acq_time"].astype(int).astype(str).str.zfill(4)
    out = pd.DataFrame({
        "latitude": d["latitude"], "longitude": d["longitude"],
        "datetime_utc": pd.to_datetime(d["acq_date"].astype(str) + " " + t, format="%Y-%m-%d %H%M", utc=True),
        "sensor": sensor, "satellite": d["satellite"].astype(str),
        "frp": d["frp"], "daynight": d["daynight"]})
    num = pd.to_numeric(d["confidence"], errors="coerce")
    if sensor == "MODIS":
        out["conf"] = num / 100
    else:
        out["conf"] = num / 100 if num.notna().all() else d["confidence"].astype(str).str.lower().map(VCONF)
    return out

def kind(name, d):
    """Decide MODIS vs VIIRS: instrument column, then VIIRS-only column, then file name, then MODIS column."""
    if "instrument" in d.columns and len(d):
        v = str(d["instrument"].iloc[0]).upper()
        if "VIIRS" in v: return "VIIRS"
        if "MODIS" in v: return "MODIS"
    if "bright_ti4" in d.columns: return "VIIRS"
    n = name.upper()
    if any(k in n for k in ("J1V", "J2V", "_SV", "VIIRS")): return "VIIRS"
    if "DL_FIRE_M" in n or "MODIS" in n: return "MODIS"
    if "bright_t31" in d.columns: return "MODIS"
    return None

def load(folder):
    mo, vi = [], []
    for name, d in read_files(folder):
        k = kind(name, d)
        if k == "MODIS": mo.append(normalize(d, "MODIS"))
        elif k == "VIIRS": vi.append(normalize(d, "VIIRS"))
        print(f"  {k or 'skipped'} : {name} ({len(d)} rows)")
    if not mo or not vi:
        sys.exit("Need at least one MODIS file and one VIIRS file in the data folder.")
    return pd.concat(mo, ignore_index=True), pd.concat(vi, ignore_index=True)

def cells(df, tag, cell):
    d = df.assign(date=df["datetime_utc"].dt.strftime("%Y-%m-%d"),
                  cx=np.floor(df["longitude"] / cell).astype(int),
                  cy=np.floor(df["latitude"] / cell).astype(int))
    return d.groupby(["date", "cx", "cy"]).agg(**{f"n_{tag}": ("frp", "size"), f"frp_{tag}": ("frp", "sum"),
                                                  f"conf_{tag}": ("conf", "max")}).reset_index()

def calendar(m, v, cell):
    c = cells(m, "modis", cell).merge(cells(v, "viirs", cell), on=["date", "cx", "cy"], how="outer")
    hm, hv = c["n_modis"].notna(), c["n_viirs"].notna()
    c["status"] = np.select([hm & hv, hm], ["both", "MODIS_only"], "VIIRS_only")
    cm, cv = c["conf_modis"].fillna(0), c["conf_viirs"].fillna(0)
    c["conf_unified"] = 1 - (1 - cm) * (1 - cv)        # simple rule: agreement raises confidence
    c["frp_total"] = c["frp_modis"].fillna(0) + c["frp_viirs"].fillna(0)
    c["lat"] = (c["cy"] + .5) * cell; c["lon"] = (c["cx"] + .5) * cell
    return c

def flagged(m, v, cell, days):
    c = calendar(m, v, cell); c = c[c["date"].isin(days)].copy()
    mr, vr = c[c.n_modis.notna()], c[c.n_viirs.notna()]
    ms, vs = set(zip(mr.date, mr.cx, mr.cy)), set(zip(vr.date, vr.cx, vr.cy))
    near = lambda s, d, x, y: any((d, x + i, y + j) in s for i in (-1, 0, 1) for j in (-1, 0, 1))
    c["viirs_nearby"] = [near(vs, d, x, y) for d, x, y in zip(c.date, c.cx, c.cy)]
    c["modis_nearby"] = [near(ms, d, x, y) for d, x, y in zip(c.date, c.cx, c.cy)]
    return c

def rates(m, v, cell, days):
    c = flagged(m, v, cell, days)
    f = lambda s: float(s.mean()) if len(s) else None
    return f(c[c.n_modis.notna()].viirs_nearby), f(c[c.n_viirs.notna()].modis_nearby)

BINS = [0, 2, 4, 8, 16, 32, 64, 1e9]; LABELS = ["<2", "2-4", "4-8", "8-16", "16-32", "32-64", "64+"]
def curve(sub, col, flag):
    g = sub.groupby(pd.cut(sub[col], BINS, labels=LABELS, right=False), observed=False)[flag].agg(["size", "mean"])
    return [{"bin": k, "n": int(r["size"]), "share": None if r["size"] == 0 else float(r["mean"])} for k, r in g.iterrows()]

def extras(m, v, cell, inner):
    c = flagged(m, v, cell, inner)
    curves = {"viirs_seen_by_modis": curve(c[c.n_viirs.notna()], "frp_viirs", "modis_nearby"),
              "modis_seen_by_viirs": curve(c[c.n_modis.notna()], "frp_modis", "viirs_nearby")}
    bench, sats = None, v["satellite"].value_counts().index.tolist()[:2]
    if len(sats) == 2:                                   # same VIIRS instrument on two satellites = best-case agreement
        A, B = v[v.satellite == sats[0]], v[v.satellite == sats[1]]
        x, y = rates(A, B, cell, inner); nx, _ = rates(A, B.assign(longitude=B["longitude"] + .5), cell, inner)
        bench = {"a": sats[0], "b": sats[1], "n_a": len(A), "n_b": len(B), "a_by_b": x, "b_by_a": y, "null": nx}
    return curves, bench

def build(m, v, cell, keep_edges):
    base = calendar(m, v, cell)
    dates = sorted(base["date"].unique())
    inner = dates if keep_edges or len(dates) <= 3 else dates[1:-1]
    for tag, flt in (("D", "D"), ("N", "N")):
        c = calendar(m[m.daynight == flt], v[v.daynight == flt], cell)[["date", "cx", "cy", "status"]]
        base = base.merge(c.rename(columns={"status": "s_" + tag}), on=["date", "cx", "cy"], how="left")
    stats = {}
    vs = v.assign(longitude=v["longitude"] + .5)           # control: shifted VIIRS = chance-level agreement
    for key, flt in (("all", None), ("day", "D"), ("night", "N")):
        mm, vv = (m, v) if flt is None else (m[m.daynight == flt], v[v.daynight == flt])
        a, b = rates(mm, vv, cell, inner)
        na, nb = rates(mm, vs if flt is None else vs[vs.daynight == flt], cell, inner)
        stats[key] = {"modis_pts": len(mm), "viirs_pts": len(vv), "modis_conf": a, "viirs_conf": b,
                      "null_modis": na, "null_viirs": nb}
    code = lambda s: s.map(CODE).fillna(0).astype(int)
    rows = [[dates.index(r.date), round(r.lat, 4), round(r.lon, 4), CODE[r.status],
             CODE.get(r.s_D, 0) if isinstance(r.s_D, str) else 0, CODE.get(r.s_N, 0) if isinstance(r.s_N, str) else 0,
             round(float(r.conf_unified), 2), round(float(r.frp_total), 1)] for r in base.itertuples()]
    curves, bench = extras(m, v, cell, inner)
    data = {"dates": dates, "inner": inner, "rows": rows, "stats": stats, "cell": cell, "curves": curves, "bench": bench,
            "bounds": [[float(base.lat.min()), float(base.lon.min())], [float(base.lat.max()), float(base.lon.max())]]}
    return base, data

def demo_files(folder):
    os.makedirs(folder, exist_ok=True); rng = np.random.default_rng(3)
    cl = [(77.8, 18.5), (79.8, 22.4), (78.2, 22.7), (77.9, 23.4), (80, 23)]
    def mk(n, viirs, sat):
        c = [cl[i] for i in rng.integers(0, len(cl), n)]
        d = pd.DataFrame({"latitude": [x[1] + rng.normal(0, .5) for x in c], "longitude": [x[0] + rng.normal(0, .5) for x in c],
            "acq_date": ["2026-11-0%d" % k for k in rng.integers(1, 8, n)],
            "acq_time": rng.integers(0, 24, n) * 100 + rng.integers(0, 60, n), "satellite": sat,
            "frp": rng.uniform(1, 40, n), "daynight": rng.choice(["D", "N"], n, p=[.8, .2])})
        if viirs: d["bright_ti4"] = 330; d["confidence"] = rng.choice(["l", "n", "h"], n)
        else: d["bright_t31"] = 300; d["confidence"] = rng.integers(0, 100, n)
        return d
    for name, d in (("DEMO_M", mk(400, False, "Terra")), ("DEMO_J1V", mk(2500, True, "N20")), ("DEMO_J2V", mk(2500, True, "N21"))):
        d.to_csv("tmp.csv", index=False)
        with zipfile.ZipFile(os.path.join(folder, name + ".zip"), "w") as z: z.write("tmp.csv", "fire.csv")
    os.remove("tmp.csv")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data"); ap.add_argument("--cell", type=float, default=0.02)
    ap.add_argument("--demo", action="store_true"); ap.add_argument("--keep-edge-days", action="store_true")
    ap.add_argument("--no-open", action="store_true"); ap.add_argument("--out", default="output")
    a = ap.parse_args()
    folder = "demo_data" if a.demo else a.data
    if a.demo: demo_files(folder); print("DEMO MODE: synthetic data, results are NOT real.")
    print("Reading", folder); m, v = load(folder)
    base, data = build(m, v, a.cell, a.keep_edge_days)
    data["demo"] = a.demo
    out = a.out; os.makedirs(out, exist_ok=True)
    base.to_csv(f"{out}/fire_calendar.csv", index=False)
    html = open("template.html", encoding="utf-8").read().replace("/*__DATA__*/null", json.dumps(data))
    open(f"{out}/index.html", "w", encoding="utf-8").write(html)
    json.dump({k: data[k] for k in ("dates", "inner", "stats", "bench", "curves", "cell")}, open(f"{out}/stats.json", "w"), indent=1)
    print(json.dumps(data["stats"], indent=1)); print("Full days used for stats:", data["inner"])
    print(f"Done. Open {out}/index.html")
    if not a.no_open: webbrowser.open("file://" + os.path.abspath(f"{out}/index.html"))

# data_prep/load_center.py
from pathlib import Path

def load_center(galaxy, *, center_root):
    path = Path(center_root) / galaxy / "center.txt"
    if not path.exists():
        raise FileNotFoundError(path)

    d = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"): continue
        k,v = line.split("=",1)
        d[k.strip()] = v.strip()

    d["ra_deg"] = float(d["ra_deg"])
    d["dec_deg"] = float(d["dec_deg"])
    d["distance_pc"] = float(d["distance_kpc"]) * 1000.0
    d["pa_deg"] = None if d["pa_deg"].lower()=="nan" else float(d["pa_deg"])

    return d

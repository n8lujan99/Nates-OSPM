from pathlib import Path
Galaxy = "Draco"
# Paths need to be specified in Data_Paths.py sticking to the following convention:
OSPM_ROOT    = Path(__file__).resolve().parent
PROFILE_ROOT = OSPM_ROOT / "data" / "profiles" / Galaxy

def build_data_paths(PROFILE_ROOT: Path, *, run_label: str = "default"):
    base = PROFILE_ROOT / run_label

    return {
        PROFILE_ROOT = OSPM_ROOT / "data" / "profiles" / Galaxy
        "DATA_CSV":       base / "data.csv",
        "STAR_CSV":       base / "star.csv",
        "OUT_CSV":        base / "ospm_input.csv",
        "PROFILE_CSV":    base / "stellar_profile.csv",
        "SPEC_PATH":      base / "spec.csv",
        "CSV_PATH":       base / "daemon_deck.csv",
        "LOG_PATH":       base / "ospm.log",
        "CHECKPOINT_DIR": base / "checkpoints",
    }


def ensure_data_dirs(paths):
    for p in paths.values():
        if isinstance(p, Path):
            if p.suffix:
                p.parent.mkdir(parents=True, exist_ok=True)
            else:
                p.mkdir(parents=True, exist_ok=True)
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
# data_prep/DATA_NEW_GAL.py
# Replaces old test_data_pipe
# Running: python -m data_prep.DATA_NEW_GAL

from OSPM.OSPM_Config import CONFIG, PROFILE_ROOT
from .Data_Paths import build_data_paths, ensure_data_dirs
from .Data_Sources import build_sources_catalog
from .Data_Preprocess import preprocess_stars_for_ospm, quality_mask_from_config, PREPROCESS_MODE
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


def main(*, run_label="default", write=True):
    paths=build_data_paths(PROFILE_ROOT,run_label=run_label); ensure_data_dirs(paths)
    df_raw=build_sources_catalog(PROFILE_ROOT=PROFILE_ROOT,scratch=not write)
    print(f"[DATA] raw sources: {len(df_raw)} rows")
    qmask=quality_mask_from_config(CONFIG) if PREPROCESS_MODE=="regular" else None

    df_clean=preprocess_stars_for_ospm(
        df_raw,
        ra_col="ra",dec_col="dec",v_col="vlos",v_err_col="vlos_err",
        ra0_deg=CONFIG["RA0_DEG"],dec0_deg=CONFIG["DEC0_DEG"],
        distance_pc=CONFIG["DISTANCE_PC"],
        r_max_pc=CONFIG.get("R_MAX_STARS_PC"),
        quality_mask=qmask,
        bin_mode=CONFIG.get("BIN_MODE","circular"),
        pa_deg=CONFIG.get("PA_DEG",0.0),
        axis_ratio_q=CONFIG.get("AXIS_RATIO_Q",1.0),
    )

    print(f"[DATA] preprocessed stars: {len(df_clean)} rows")
    if write:
        df_clean.to_csv(paths["STAR_CSV"],index=False)
        print(f"[OK] wrote STAR_CSV -> {paths['STAR_CSV']}")
    return df_clean
if __name__=="__main__": main(run_label="default",write=True)



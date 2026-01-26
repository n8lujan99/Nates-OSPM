"""
OSPM observational preprocessing.
Pure geometry + irreversible cuts.
"""
import numpy as np
from .Data_Geometry import project_tangent_plane, elliptical_coordinates, circular_radius
PREPROCESS_MODE="dwarf"
def preprocess_defaults(mode):
    if mode=="dwarf":   return dict(use_quality_mask=False,require_v_err=False,r_max_pc=None)
    if mode=="regular": return dict(use_quality_mask=True, require_v_err=True, r_max_pc=None)
    raise ValueError(f"Unknown preprocess mode: {mode}")
def quality_mask_from_config(cfg):
    def mask(df):
        m=np.ones(len(df),bool)
        if "ruwe" in df and "RUWE_MAX" in cfg: m&=(df["ruwe"]<=cfg["RUWE_MAX"])
        if "parallax" in df and "parallax_error" in df and "PAR_SNR_MIN" in cfg:
            snr=df["parallax"]/df["parallax_error"]
            m&=((snr>=cfg["PAR_SNR_MIN"])|(~np.isfinite(snr)))
        return m
    return mask

def preprocess_stars_for_ospm(
    df, *,
    ra_col="ra", dec_col="dec",
    v_col="vlos", v_err_col="vlos_err",
    ra0_deg=None, dec0_deg=None, distance_pc=None,
    x_col=None, y_col=None,
    r_max_pc=None, quality_mask=None,
    bin_mode="circular",
    pa_deg=0.0, axis_ratio_q=1.0,
    mode=PREPROCESS_MODE
):
    d = preprocess_defaults(mode)
    if not d["use_quality_mask"]:
        quality_mask = None
    if not d["require_v_err"]:
        v_err_col = None
    if r_max_pc is None:
        r_max_pc = d["r_max_pc"]

    df = df.copy()

    # --------------------------------------------------
    # Geometry
    # --------------------------------------------------
    if x_col is None or y_col is None:
        if ra0_deg is None or dec0_deg is None or distance_pc is None:
            raise ValueError("ra0_deg, dec0_deg, and distance_pc required")

        x_pc, y_pc, r_pc = project_tangent_plane(
            df[ra_col].to_numpy(),
            df[dec_col].to_numpy(),
            ra0_deg, dec0_deg,
            distance_pc
        )
        df["x_pc"] = x_pc
        df["y_pc"] = y_pc
        df["r_pc"] = r_pc
    else:
        df["x_pc"] = df[x_col].to_numpy()
        df["y_pc"] = df[y_col].to_numpy()
        df["r_pc"] = circular_radius(df["x_pc"], df["y_pc"])

    # --------------------------------------------------
    # Radial coordinate
    # --------------------------------------------------
    if bin_mode == "elliptical":
        R, _ = elliptical_coordinates(
            df["x_pc"], df["y_pc"],
            pa_deg=pa_deg,
            axis_ratio_q=axis_ratio_q
        )
        df["R_ell"] = R
    else:
        R = df["r_pc"].to_numpy()

    # --------------------------------------------------
    # Velocities
    # --------------------------------------------------
    df[v_col] = df[v_col].replace(-9.999, np.nan)
    if v_err_col is not None and v_err_col in df:
        df[v_err_col] = df[v_err_col].replace(-9.999, np.nan)

    df["has_vlos"] = np.isfinite(df[v_col]) & (
        np.isfinite(df[v_err_col])
        if v_err_col is not None and v_err_col in df
        else True
    )

    # --------------------------------------------------
    # Masks
    # --------------------------------------------------
    m = np.isfinite(R)
    if quality_mask is not None:
        m &= quality_mask(df)
    if r_max_pc is not None:
        m &= (R <= r_max_pc)

    out = df.loc[m].reset_index(drop=True)
    out["preprocess_mode"] = mode
    return out

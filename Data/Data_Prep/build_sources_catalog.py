import pandas as pd
import numpy as np

def load_star_catalog(path, *, distance_pc):
    cols = [
        "id",
        "vlos",
        "vlos_err",
        "ew",
        "ew_err",
        "r_arcmin",
        "mjd",
        "g_mag",
        "r_mag",
        "i_mag",
        "member_flag",
        "p_em",
        "p_bayes",
    ]

    df = pd.read_csv(
        path,
        delim_whitespace=True,
        names=cols,
        comment="#"
    )

    # Sentinel cleanup
    for c in ["vlos", "vlos_err", "p_em", "p_bayes"]:
        df[c] = df[c].replace([-9.999, 0.0], np.nan)

    # Projected radius → physical radius
    df["r_pc"] = distance_pc * np.deg2rad(df["r_arcmin"] / 60.0)

    # OSPM-required flags
    df["has_vlos"] = np.isfinite(df["vlos"]) & np.isfinite(df["vlos_err"])

    # This is a radial catalog, not a sky catalog
    df["ra"]  = np.nan
    df["dec"] = np.nan
    
    print("rows:", len(df))
    print("finite vlos:", np.isfinite(df["vlos"]).sum())
    print("finite vlos_err:", np.isfinite(df["vlos_err"]).sum())
    print(df.loc[np.isfinite(df["vlos"]), ["id","vlos","vlos_err","r_arcmin"]].head(5))

    return df


import os
import pandas as pd
import numpy as np

def load_pipeline_state(repo_root="Nates-OSPM"):

    data_root = os.path.join(repo_root, "Data")

    # Find first Gal_* directory
    gal_dirs = [d for d in os.listdir(data_root)
                if d.startswith("Gal_") and
                os.path.isdir(os.path.join(data_root, d))]

    if not gal_dirs:
        raise RuntimeError("No Gal_* directories found.")

    gal_path = os.path.join(data_root, gal_dirs[0],
                            "whichgalaxy", "default")

    data_file = os.path.join(gal_path, "data.csv")
    deck_file = os.path.join(gal_path, "daemon_deck.csv")

    if not os.path.exists(data_file):
        raise RuntimeError("data.csv not found.")
    if not os.path.exists(deck_file):
        raise RuntimeError("daemon_deck.csv not found.")

    # Load star data
    df_stars = pd.read_csv(data_file)

    stars = {
        "R_proj": df_stars["R_star_m"].values,
        "v_los": df_stars["v_star_mps"].values,
        "has_vlos": df_stars["has_vlos"].astype(bool).values
    }

    # Load daemon deck
    df_deck = pd.read_csv(deck_file)

    # Column 4 = chi2 (zero-indexed 3)
    chi_col = df_deck.columns[3]

    best_idx = df_deck[chi_col].idxmin()
    best_row = df_deck.loc[best_idx]

    # Assuming first 3 columns are theta
    rho_s = best_row[0]
    r_s   = best_row[1]
    MBH   = best_row[2]

    theta_best = (rho_s, r_s, MBH)

    print("Best-fit theta:", theta_best)
    print("Minimum chi²:", best_row[chi_col])

    return theta_best, stars

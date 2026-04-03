#!/bin/bash

cd /home/nate/research/Nates-OSPM
cat <<'EOF' >> Data/Gal_Profiles/Draco/Draco_OSPM_Config.py
# --- LOM suggested region update for Draco ---
# generated: 2026-03-29T14:42:47
"PARAMETER_NAMES": ["rho_s", "r_s", "MBH"],
"INITIAL_THETA": [6.327, 3436.0, 2.250e+06],
"THETA_BOUNDS": [
   (0.6013, 14.82),      # rho_s
   (50.0, 7987.0),    # r_s
   (3554.0, 5.000e+06),      # MBH
],
EOF

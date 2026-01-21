# Data_Geometry.py
# PURE GEOMETRY UTILITIES FOR OSPM PREPROCESSING
# Should not need frequent edits
"""
Geometric transforms for OSPM preprocessing.
Pure geometry. No statistics. No binning. No IO.
"""
import numpy as np

def project_tangent_plane(ra_deg, dec_deg, ra0_deg, dec0_deg, distance_pc):
    ra = np.deg2rad(ra_deg); dec = np.deg2rad(dec_deg)
    ra0 = np.deg2rad(ra0_deg); dec0 = np.deg2rad(dec0_deg)
    dra = (ra - ra0) * np.cos(dec0); ddec = dec - dec0
    x_pc = distance_pc * dra; y_pc = distance_pc * ddec
    r_pc = np.hypot(x_pc, y_pc)
    return x_pc, y_pc, r_pc

def elliptical_coordinates(x_pc, y_pc, pa_deg=0.0, axis_ratio_q=1.0):
    pa = np.deg2rad(pa_deg); c = np.cos(pa); s = np.sin(pa)
    x = np.asarray(x_pc, float); y = np.asarray(y_pc, float)
    x_maj = c * x + s * y; y_min = -s * x + c * y
    R_ell = np.sqrt(x_maj**2 + (y_min / axis_ratio_q)**2)
    theta = np.arctan2(y_min / axis_ratio_q, x_maj)
    return R_ell, theta

def circular_radius(x_pc, y_pc):
    x = np.asarray(x_pc, float); y = np.asarray(y_pc, float)
    return np.hypot(x, y)

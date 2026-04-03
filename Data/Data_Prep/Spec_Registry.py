# Data/Data_Prep/Spec_Registry.py
# Known Vizier catalogs with stellar radial velocities for dwarf galaxies.
# Each entry: catalog_id, galaxy filter (None = all rows belong to this galaxy),
# column mappings for ra/dec/vlos/vlos_err, and an optional Target column + pattern.

# Format per entry:
#   {
#       "catalog":   "J/ApJS/268/19/...",
#       "table":     "tablename" or None (first table),
#       "ra_col":    column name for RA  (degrees),
#       "dec_col":   column name for Dec (degrees),
#       "v_col":     column name for radial velocity,
#       "v_err_col": column name for velocity error (or None),
#       "target_col": column to filter on (or None = use all rows),
#       "target_pat": regex pattern to match in target_col (or None),
#       "ref":       short citation for provenance
#   }

# ---- mega-catalog: Walker+2023 (M2FS + Hectochelle) ----
def _walker2023(galaxy_pat, table="hectocat"):
    tbl = {"hectocat": "J/ApJS/268/19/hectocat",
           "m2fshi": "J/ApJS/268/19/m2fshi",
           "m2fsmed": "J/ApJS/268/19/m2fsmed"}
    return {
        "catalog": tbl[table],
        "ra_col": "RAJ2000", "dec_col": "DEJ2000",
        "v_col": "VlosAvg", "v_err_col": "e_VlosAvg",
        "target_col": "Target", "target_pat": galaxy_pat,
        "ref": f"Walker+2023 ({table})",
    }

# ---- Walker+2009 (Carina, Fornax, Sculptor, Sextans) ----
def _walker2009(galaxy_pat):
    return {
        "catalog": "J/AJ/137/3100/stars",
        "ra_col": "RAJ2000", "dec_col": "DEJ2000",
        "v_col": "<HV>", "v_err_col": "e_<HV>",
        "target_col": "Target", "target_pat": galaxy_pat,
        "ref": "Walker+2009",
    }

# ============================================================
# REGISTRY: galaxy name -> list of catalog specs
# ============================================================
SPEC_REGISTRY = {

    # --- Bootes I ---
    "BootesI": [
        _walker2023(r"(?i)boo", "hectocat"),
        {
            "catalog": "J/ApJ/736/146/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RV", "v_err_col": "e_RV",
            "target_col": None, "target_pat": None,
            "ref": "Koposov+2011",
        },
        {
            "catalog": "J/ApJ/920/92/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RVel", "v_err_col": "e_RVel",
            "target_col": None, "target_pat": None,
            "ref": "Jenkins+2021",
        },
    ],

    # --- Carina ---
    "Carina": [
        _walker2023(r"(?i)carina|car", "m2fshi"),
        _walker2009(r"(?i)^Car"),
        {
            "catalog": "J/ApJ/830/126/table2",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RVel", "v_err_col": "e_RVel",
            "target_col": None, "target_pat": None,
            "ref": "Fabrizio+2016",
        },
    ],

    # --- Carina II ---
    "CarinaII": [
        {
            "catalog": "J/ApJ/857/145/table2",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Li+2018",
        },
    ],

    # --- Carina III ---
    "CarinaIII": [
        {
            "catalog": "J/ApJ/857/145/table3",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Li+2018",
        },
    ],

    # --- CVnI (Canes Venatici I) ---
    "CVnI": [
        _walker2023(r"(?i)cvn|canes", "hectocat"),
    ],

    # --- Coma Berenices ---
    "ComBer": [
        {
            "catalog": "J/ApJ/708/560/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RVH", "v_err_col": None,
            "target_col": "Gal", "target_pat": r"(?i)com",
            "ref": "Frebel+2010",
        },
    ],

    # --- Draco ---
    "Draco": [
        _walker2023(r"(?i)draco|dra", "hectocat"),
        {
            "catalog": "J/AJ/156/257/table2m",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "vlosmean", "v_err_col": "e_vlosmean",
            "target_col": None, "target_pat": None,
            "ref": "Spencer+2018",
        },
    ],

    # --- Draco II ---
    "DracoII": [
        {
            "catalog": "J/MNRAS/458/L59/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RV", "v_err_col": "e_RV",
            "target_col": None, "target_pat": None,
            "ref": "Martin+2016",
        },
    ],

    # --- Fornax ---
    "Fornax": [
        _walker2023(r"(?i)fornax|for", "m2fshi"),
        _walker2009(r"(?i)^For"),
        {
            "catalog": "J/ApJ/923/77/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "VLOS", "v_err_col": "e_VLOS",
            "target_col": None, "target_pat": None,
            "ref": "Pace+2021",
        },
        {
            "catalog": "J/AJ/131/2114/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Walker+2006",
        },
    ],

    # --- Grus I ---
    "GrusI": [
        _walker2023(r"(?i)grus|gru", "m2fshi"),
    ],

    # --- Grus II ---
    "GrusII": [
        {
            "catalog": "J/ApJ/892/137/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": "Sat", "target_pat": r"(?i)gru.*II|gru.*2",
            "ref": "Simon+2020",
        },
    ],

    # --- Hercules ---
    "Hercules": [
        # Walker+2023 may have some; add if confirmed
    ],

    # --- Hydrus I ---
    "HydrusI": [
        _walker2023(r"(?i)hydrus|hyi", "m2fshi"),
    ],

    # --- Leo I ---
    "LeoI": [
        {
            "catalog": "J/ApJ/675/201/table5",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HV", "v_err_col": None,
            "target_col": None, "target_pat": None,
            "ref": "Mateo+2008",
        },
        {
            "catalog": "J/ApJ/657/241/table2",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Koch+2007",
        },
    ],

    # --- Leo II ---
    "LeoII": [
        {
            "catalog": "J/AJ/153/254/table3",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RV", "v_err_col": "e_RV",
            "target_col": None, "target_pat": None,
            "ref": "Spencer+2017",
        },
    ],

    # --- Leo IV ---
    "LeoIV": [
        {
            "catalog": "J/ApJ/920/92/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RVel", "v_err_col": "e_RVel",
            "target_col": None, "target_pat": None,
            "ref": "Jenkins+2021",
        },
    ],

    # --- Leo V ---
    "LeoV": [
        {
            "catalog": "J/ApJ/885/53/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "Vel", "v_err_col": "e_Vel",
            "target_col": None, "target_pat": None,
            "ref": "Mutlu-Pakdil+2019",
        },
    ],

    # --- Reticulum II ---
    "ReticulumII": [
        _walker2023(r"(?i)ret", "m2fshi"),
        {
            "catalog": "J/ApJ/808/95/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "Vel", "v_err_col": "e_Vel",
            "target_col": None, "target_pat": None,
            "ref": "Simon+2015",
        },
    ],

    # --- Sculptor ---
    "Sculptor": [
        _walker2023(r"(?i)sculptor|scl", "m2fshi"),
        _walker2009(r"(?i)^Scl"),
    ],

    # --- Segue1 ---
    "Segue1": [
        _walker2023(r"(?i)seg.*1|segue", "hectocat"),
        {
            "catalog": "J/ApJ/733/46/table3",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "Vel", "v_err_col": "e_Vel",
            "target_col": None, "target_pat": None,
            "ref": "Simon+2011",
        },
        {
            "catalog": "J/ApJ/692/1464/table2",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Geha+2009",
        },
    ],

    # --- Segue 2 ---
    "Segue2": [
        {
            "catalog": "J/ApJ/770/16/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Kirby+2013",
        },
    ],

    # --- Sextans ---
    "Sextans": [
        _walker2023(r"(?i)sextans|sex", "hectocat"),
        _walker2023(r"(?i)sextans|sex", "m2fshi"),
        _walker2009(r"(?i)^Sex"),
    ],

    # --- Tucana II ---
    "TucanaII": [
        _walker2023(r"(?i)tuc.*2|tuc.*II", "m2fshi"),
        {
            "catalog": "J/ApJ/838/83/table2",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Simon+2019 (Walker 2016 reanalysis)",
        },
    ],

    # --- Tucana III ---
    "TucanaIII": [
        {
            "catalog": "J/ApJ/838/11/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": None, "target_pat": None,
            "ref": "Simon+2017",
        },
    ],

    # --- Tucana IV, Tucana V, Grus II (Simon+2020) ---
    "TucanaIV": [
        {
            "catalog": "J/ApJ/892/137/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "HRV", "v_err_col": "e_HRV",
            "target_col": "Sat", "target_pat": r"(?i)tuc.*IV|tuc.*4",
            "ref": "Simon+2020",
        },
    ],

    # --- Ursa Major II ---
    "UMaII": [
        _walker2023(r"(?i)uma.*II|uma.*2|ursa.*maj", "hectocat"),
        {
            "catalog": "J/ApJ/708/560/table1",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RVH", "v_err_col": None,
            "target_col": "Gal", "target_pat": r"(?i)uma",
            "ref": "Frebel+2010",
        },
    ],

    # --- Ursa Minor ---
    "UMi": [
        _walker2023(r"(?i)umi|ursa.*min", "hectocat"),
        {
            "catalog": "J/AJ/156/257/table4",
            "ra_col": "RAJ2000", "dec_col": "DEJ2000",
            "v_col": "RV", "v_err_col": "e_RV",
            "target_col": None, "target_pat": None,
            "ref": "Spencer+2018",
        },
    ],

    # --- Crater 2 ---
    "Crater2": [
        _walker2023(r"(?i)crater|crt", "hectocat"),
    ],

    # --- Phoenix 2 ---
    "Phoenix2": [
        _walker2023(r"(?i)phoe|phx", "m2fshi"),
    ],

    # --- Antlia 2 ---
    "Antlia2": [
        _walker2023(r"(?i)antlia|ant", "m2fshi"),
    ],

    # --- Crater 1 ---
    "Crater1": [
        _walker2023(r"(?i)crater.*1|crt.*1", "m2fshi"),
    ],
}

# OSPM Data Pipeline — README

This directory defines the **only valid observational data path into OSPM**.

Its job is to take heterogeneous astronomical catalogs, turn them into a single frozen observational snapshot, and then reduce that snapshot to a galaxy-anchored stellar sample suitable for dynamical modeling.

Nothing more. Nothing less.

---

## What this pipeline does

### 1. Harvests stars from the sky, not from models
The pipeline queries multiple external catalogs (Gaia, SDSS, LAMOST, APOGEE, etc.) using a **declared sky position and radius**.

Each source is queried independently.  
No source is trusted more than another at ingest time.  
No catalog decides membership.

The only shared assumption at this stage is geometry:
same center, same angular radius.

---

### 2. Fuses measurements without inventing information
Stars from different catalogs are matched by sky position.

For each matched star:
- velocities are compared
- the measurement with the **smallest fractional error** is selected
- the source of that measurement is recorded

No averaging.  
No fitting.  
No inference.

If a star has no velocity, it does not survive fusion.

---

### 3. Writes a frozen observational CSV
After fusion, the pipeline writes a CSV file that represents:

the maximum verified observational state of knowledge for this galaxy under this configuration.

This CSV is treated as **immutable**.

Downstream code may read it.  
Downstream code may transform it in memory.  
Downstream code must never modify it in place.

If something is wrong, the pipeline is rerun.

---

### 4. Anchors stars in a declared galaxy frame
In preprocessing, stars are projected into a flat galaxy-centric coordinate system using:

- a declared RA/Dec center  
- a declared distance  

These values come from configuration or literature, not from the data.

Velocities are recentred by subtracting the systemic median.
This places the galaxy at rest by construction.

---

### 5. Applies hard, irreversible cuts
Stars are removed if they fail basic, non-negotiable criteria:

- missing or non-finite positions or velocities  
- missing velocity errors (if required)  
- outside the declared maximum stellar radius  
- optional quality cuts (RUWE, parallax SNR, etc.)

Once a star is removed here, it is gone forever.

No probabilistic membership.  
No later rescue.

---

## What this pipeline does not do

### It does not infer galaxy properties
The pipeline does not:
- fit the galaxy center
- infer ellipticity
- infer scale radii
- infer membership probabilities

Galaxy geometry is **declared**, not learned.

This avoids feedback loops where:
- foreground stars move the center
- membership changes the coordinate system
- the solver chases a moving target

---

### It does not bin, fit, or model
There is:
- no binning
- no dispersion calculation
- no likelihood evaluation
- no orbit integration
- no χ²

Those belong to OSPM proper.

This pipeline stops at “clean stars in a galaxy frame”.

---

### It does not mutate data downstream
Once the frozen CSV is written:
- it is never edited
- it is never appended
- it is never partially rewritten

All downstream work operates on in-memory copies.

This guarantees reproducibility.

---

## Why the design looks the way it does

### Why multiple source functions with no arguments
Each source function is a **harvester**, not a transformer.

Its real inputs are:
- the declared sky geometry
- the external catalog itself

Passing DataFrames into source functions would bias ingestion and collapse pipeline stages.

Harvest first.  
Fuse second.

---

### Why the center is declared, not computed
Computing the center from stars couples geometry to membership.

That breaks dynamical modeling.

The center may be diagnosed from data, but OSPM always obeys a declared center.

---

### Why quality masks are optional
Quality criteria are **policy**, not physics.

The pipeline allows:
- no quality mask
- a Gaia-only mask
- a custom survey-specific mask

The preprocess function does not assume one exists.

---

### Why this pipeline is intentionally conservative
OSPM is sensitive to:
- geometry
- boundaries
- consistency

It is safer to throw away borderline stars early than to explain bad fits later.

This pipeline prefers:
- fewer stars
- cleaner assumptions
- explicit decisions

---

## Mental model to keep in mind

Think of the pipeline as drawing a line in the sand.

Everything before the line:
catalogs, networks, mess, ambiguity.

Everything after the line:
fixed data, fixed geometry, fixed stars, deterministic modeling.

If you ever feel tempted to “just tweak the stars later”, you are on the wrong side of that line.

---

## One-sentence summary

This pipeline freezes the observational universe so the dynamical solver can stop worrying about it.

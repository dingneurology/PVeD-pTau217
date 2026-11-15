"""Run the full Python analysis pipeline in order.

    python run_all.py            # everything
    python run_all.py --from prep.py   # resume from a given step

Then draw Figures 1 (data panels), 2, 3 and Supplementary Figure S1 in R:
    Rscript figures/R/make_figures.R
"""
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
STEPS = [
    ("analysis/build_adni.py", "ADNI: merge baseline file with PET, MRI and CDR closest to the PVeD MRI"),
    ("analysis/build_habs.py", "HABS-HD: merge PVeD, plasma, PET, MRI, clinical and genomic releases"),
    ("analysis/prep.py", "Harmonised analytic datasets, standardisation and pTau217-PVeD profiles"),
    ("analysis/analyses.py", "Cross-sectional models, sensitivity analyses, Table 1, longitudinal mixed models"),
    ("analysis/lmm_profile.py", "Profile-by-time mixed models of PACC change"),
    ("analysis/incr.py", "Incremental explained variance"),
    ("analysis/extra.py", "Additional sensitivity analyses"),
    ("analysis/tau_tracer.py", "ADNI tau PET tracer handling"),
    ("analysis/long_methods.py", "ADNI longitudinal model-form comparison"),
    ("analysis/meta_long.py", "Two-cohort random-effects synthesis of longitudinal estimates"),
    ("figures/figs.py", "Figure 4 and Supplementary Figures S2-S3"),
    ("figures/export_fig_sources.py", "Source data for the R figures"),
]

if __name__ == "__main__":
    start = 0
    if "--from" in sys.argv:
        name = sys.argv[sys.argv.index("--from") + 1]
        start = next(i for i, (s, _) in enumerate(STEPS) if s.endswith(name))
    for script, what in STEPS[start:]:
        t0 = time.time()
        print(f"\n=== {script}: {what}", flush=True)
        subprocess.run([sys.executable, os.path.join(ROOT, script)], check=True, cwd=ROOT)
        print(f"--- done in {time.time() - t0:.0f} s", flush=True)
    print("\nAll Python steps finished. Results: results/. Now run: Rscript figures/R/make_figures.R")

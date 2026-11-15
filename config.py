"""Central path configuration.

Edit the three input locations below, or set the environment variables
PVED_ADNI_DIR, PVED_HABS_DIR and PVED_PHC_TAU_FILE. All other folders are created inside the repository.
"""
import os

ROOT = os.path.dirname(os.path.abspath(__file__))


def _dir(env, default):
    return os.path.join(os.path.abspath(os.environ.get(env, default)), "")


# ---- inputs (not distributed; obtain from ADNI / HABS-HD under their data use agreements) ----
ADNI = _dir("PVED_ADNI_DIR", os.path.join(ROOT, "data", "ADNI_PVED"))
HABS = _dir("PVED_HABS_DIR", os.path.join(ROOT, "data", "HABS_PVED"))
PHC_TAU = os.path.abspath(os.environ.get("PVED_PHC_TAU_FILE", os.path.join(ROOT, "data", "ADSP_PHC_PET_Tau_Detailed.csv")))

# ---- generated ----
WORK = _dir("PVED_WORK_DIR", os.path.join(ROOT, "work"))                 # intermediate participant-level files
OUT = _dir("PVED_OUT_DIR", os.path.join(ROOT, "results"))                # result tables (CSV)
FIG = os.path.join(ROOT, "figures", "output", "")                         # Python figures (Figure 4, Supplementary S2-S3)
FIG_SOURCE = os.path.join(ROOT, "figures", "R", "source_data", "")        # inputs for figures/R/make_figures.R

for _p in (WORK, OUT, FIG, FIG_SOURCE):
    os.makedirs(_p, exist_ok=True)

from pkg_resources import get_distribution

try:
    __version__ = get_distribution("varifier").version
except:
    __version__ = "local"


__all__ = [
    "dnadiff",
    "edit_distance",
    "probe",
    "probe_mapping",
    "recall",
    "tasks",
    "truth_variant_finding",
    "utils",
    "vcf_evaluate",
    "vcf_stats",
]

from varifier import *

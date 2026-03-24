from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
INPUT_TSV = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"
MODELS_DIR = PROJECT_ROOT / "cif_models"


# Optional: set this to a UniProt ID like "O00571" to limit step 2.
RUN_ONLY_UNIPROT: str | None = None


def run_step(cmd: list[str]) -> None:
    print("\n$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    print("=== Bio465 Capstone Pipeline ===")
    python_exe = sys.executable

    step1 = [python_exe, str(SCRIPTS_DIR / "1_filter.py")]

    step2 = [
        python_exe,
        str(SCRIPTS_DIR / "2_download_structures.py"),
        str(INPUT_TSV),
        "--id_column",
        "uniprot_id",
        "--out_dir",
        str(MODELS_DIR),
        "--prefer",
        "cif",
        "--delay",
        "0.1",
        "--also_pae",
    ]

    step3 = [python_exe, str(SCRIPTS_DIR / "3_find_nearby_mutations.py")]
    if RUN_ONLY_UNIPROT:
        step3.extend(["--uniprot", RUN_ONLY_UNIPROT])

    print("Running pipeline steps in order:")
    print("1) Filter and merge PTMD + TCGA data")
    print("2) Download AlphaFold CIF models and PAE files")
    print("3) Find nearby mutations and compute distances")

    run_step(step1)
    run_step(step2)
    run_step(step3)

    print("\nPipeline complete.")


if __name__ == "__main__":
	main()

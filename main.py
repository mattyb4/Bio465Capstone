from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
INPUT_TSV = PROJECT_ROOT / "data" / "PTMD_TCGA_hotspots_by_protein.tsv"
MODELS_DIR = PROJECT_ROOT / "cif_models"


# Optional: set this to a UniProt ID like "O00571" to limit step 2.
RUN_ONLY_UNIPROT: str | None = None


def run_step(cmd: list[str]) -> None:
    print("\n$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    python_exe = sys.executable

    step1 = [
        python_exe,
        str(SCRIPTS_DIR / "getAlphafoldCifs.py"),
        str(INPUT_TSV),
        "--id_column",
        "uniprot_id",
        "--out_dir",
        str(MODELS_DIR),
        "--prefer",
        "cif",
        "--delay",
        "0.1",
    ]

    step2 = [python_exe, str(SCRIPTS_DIR / "create_nearby_mutations_db.py")]
    if RUN_ONLY_UNIPROT:
        step2.extend(["--uniprot", RUN_ONLY_UNIPROT])

    step3 = [python_exe, str(SCRIPTS_DIR / "calcDist.py")]

    print("Running pipeline steps in order:")
    print("1) Download CIF models")
    print("2) Build nearby_mutations_db.tsv")
    print("3) Build ptm_linear_distances.tsv")

    run_step(step1)
    run_step(step2)
    run_step(step3)

    print("\nPipeline complete.")


if __name__ == "__main__":
	main()

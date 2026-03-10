import numpy as np
from pathlib import Path
import re
import csv
import argparse
from typing import Any
from biotite.structure.io.pdbx import CIFFile, get_structure  # type: ignore[import-untyped]

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
MODELS_ROOT = PROJECT_ROOT / "cif_models"
OUTPUT_PATH = PROJECT_ROOT / "Output" / "nearby_mutations_db.tsv"
PTM_TSV_PATH = PROJECT_ROOT / "data" / "PTMD_TCGA_hotspots_by_protein.tsv"

_PTM_ROWS: list[dict[str, Any]] | None = None


def get_ptm_rows():
    global _PTM_ROWS
    if _PTM_ROWS is None:
        with PTM_TSV_PATH.open("r", encoding="utf-8", newline="") as handle:
            _PTM_ROWS = list(csv.DictReader(handle, delimiter="\t"))
    return _PTM_ROWS

def get_ca_coord(chain, residue_number):
    mask = (chain.res_id == residue_number) & (chain.atom_name == "CA")
    if not np.any(mask):
        return None
    return chain.coord[mask][0]

#compute distance between two 3D points
def compute_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

#find mutations within cutoff distance of PTM site
def find_nearby_mutations(chain, ptm_pos, mutation_entries, cutoff=10.0): #adjust cutoff as needed
    results = []

    ptm_coord = get_ca_coord(chain, ptm_pos)

    if ptm_coord is None:
        print(f"PTM residue {ptm_pos} not found in structure.")
        return results

    for mutation, mut_pos in mutation_entries:
        mut_coord = get_ca_coord(chain, mut_pos)

        if mut_coord is None:
            continue

        distance = compute_distance(ptm_coord, mut_coord)

        if distance <= cutoff:
            results.append({
                "mutation": mutation,
                "mutation_pos": mut_pos,
                "distance": distance
            })

    return results


#extracts PTM positions/types from input db
PTM_RE = re.compile(r"([A-Z])(\d+)")  # e.g., S557


def parse_ptm_entries(uniprot):
    entries = set()

    for row in get_ptm_rows():
        if row.get("uniprot_id") != uniprot:
            continue
        field = row.get("ptms_on_protein", "")
        for token in re.split(r"[;,]", field):
            token = token.strip()
            match = PTM_RE.search(token)
            if match:
                ptm_type = f"{match.group(1)}{match.group(2)}"
                entries.add((ptm_type, int(match.group(2))))

    return sorted(entries, key=lambda x: x[1])


def parse_gene_name(uniprot):
    for row in get_ptm_rows():
        if row.get("uniprot_id") == uniprot:
            return row.get("gene", "")

    return ""

MUT_RE = re.compile(r"([A-Z])(\d+)([A-Z*])")  # e.g., R482H
#extracts mutation positions from input db and puts them into a list
def parse_mutation_positions(uniprot=None):
    mutation_entries = set()

    for row in get_ptm_rows():
        if uniprot and row.get("uniprot_id") != uniprot:
            continue

        fields = [row.get("mutations_on_protein", "")]
        for field in fields:
            for token in re.split(r"[;,]", field):
                match = MUT_RE.search(token.strip())
                if match:
                    mutation = f"{match.group(1)}{match.group(2)}{match.group(3)}"
                    mutation_entries.add((mutation, int(match.group(2))))

    return sorted(mutation_entries, key=lambda x: (x[1], x[0]))


def find_model_file(uniprot_dir):
    candidates = sorted(uniprot_dir.glob("*.cif"))
    return candidates[0] if candidates else None


def load_first_chain(model_file):
    try:
        cif = CIFFile.read(str(model_file))
        structure = get_structure(cif, model=1)
    except Exception as exc:
        print(f"Skipping {model_file.name}: failed to parse ({exc}).")
        return None

    if structure is None or len(structure) == 0:
        print(f"Skipping {model_file.name}: no models found.")
        return None

    chain_ids = list(dict.fromkeys(structure.chain_id))
    if not chain_ids:
        print(f"Skipping {model_file.name}: no chains found.")
        return None

    chain_id = chain_ids[0]
    return structure[structure.chain_id == chain_id]


def format_mutations(hits):
    if not hits:
        return ""
    parts = [f"{hit['mutation']}-{hit['distance']:.2f}Å" for hit in sorted(hits, key=lambda h: (h["mutation_pos"], h["mutation"]))]
    return ", ".join(parts)

#This is for debugging specific cases. Run with --uniprot P12345 to only process that UniProt ID.
#Can probably be removed for final version
parser = argparse.ArgumentParser(description="Scan AFDB models for nearby PTM mutations.")
parser.add_argument("--uniprot", help="Limit processing to a single UniProt ID.")
args = parser.parse_args()


# Main processing loop: iterate over AlphaFold models, find nearby mutations for each PTM site, and write results to TSV.
with OUTPUT_PATH.open("w", encoding="utf-16", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "UniProt",
        "gene",
        "ptm_type",
        "ptm_pos",
        "mutations_within_5_positions",
        "mutation_count_within_5_positions",
        "mutations_more_than_5_positions",
        "mutation_count_more_than_5_positions",
    ])

    for uniprot_dir in sorted(MODELS_ROOT.iterdir()):
        if not uniprot_dir.is_dir():
            continue

        uniprot = uniprot_dir.name
        if args.uniprot and uniprot != args.uniprot:
            continue
        gene = parse_gene_name(uniprot)
        ptm_entries = parse_ptm_entries(uniprot)
        if not ptm_entries:
            continue
        mutation_entries = parse_mutation_positions(uniprot=uniprot)
        if not mutation_entries:
            continue

        model_file = find_model_file(uniprot_dir)
        if model_file is None:
            print(f"Skipping {uniprot}: no CIF file found.")
            continue

        chain = load_first_chain(model_file)
        if chain is None:
            continue

        for ptm_type, ptm_position in ptm_entries:
            nearby = find_nearby_mutations(chain, ptm_position, mutation_entries)
            if not nearby:
                continue
            within_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) <= 5]
            beyond_5 = [hit for hit in nearby if abs(hit["mutation_pos"] - ptm_position) > 5]
            writer.writerow([
                uniprot,
                gene,
                ptm_type,
                ptm_position,
                format_mutations(within_5),
                len(within_5),
                format_mutations(beyond_5),
                len(beyond_5),
            ])

print(f"Wrote nearby mutation data to {OUTPUT_PATH}")

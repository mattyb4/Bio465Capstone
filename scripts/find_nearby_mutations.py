from Bio.PDB import MMCIFParser
import numpy as np
from pathlib import Path
import re
import csv

parser = MMCIFParser(QUIET=True)
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
model_file = PROJECT_ROOT / "AlphaFold_models" / "P35222" / "P35222_model_0.cif" #change to be path to desired .cif file
structure = parser.get_structure("protein", model_file)

model = structure[0]   # First (and only) model
chain = list(model.get_chains())[0]  # AlphaFold usually has one chain

def get_ca_coord(chain, residue_number):
    """
    Returns CA coordinates for a given residue number.
    """
    for residue in chain:
        if residue.get_id()[1] == residue_number:
            if "CA" in residue:
                return residue["CA"].get_coord()
    return None

#compute distance between two 3D points
def compute_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

#find mutations within cutoff distance of PTM site
def find_nearby_mutations(chain, ptm_pos, mutation_positions, cutoff=10.0): #adjust cutoff as needed
    results = []

    ptm_coord = get_ca_coord(chain, ptm_pos)

    if ptm_coord is None:
        print(f"PTM residue {ptm_pos} not found in structure.")
        return results

    for mut_pos in mutation_positions:
        mut_coord = get_ca_coord(chain, mut_pos)

        if mut_coord is None:
            continue

        distance = compute_distance(ptm_coord, mut_coord)

        if distance <= cutoff:
            results.append({
                "mutation_pos": mut_pos,
                "distance": distance
            })

    return results

MUT_RE = re.compile(r"[A-Z](\d+)[A-Z*]")  # e.g., R482H, S2054L

def parse_ptm_positions(uniprot):
    tsv_path = PROJECT_ROOT / "data" / "PTM_Associated_By_PTM_PrevalenceFilteredFar.tsv"
    positions = set()

    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("UniProt") == uniprot:
                positions.add(int(row["ptm_pos"]))

    return sorted(positions)

def parse_mutation_positions(ptm_pos, uniprot=None): 
    tsv_path = PROJECT_ROOT / "data" / "PTM_Associated_By_PTM_PrevalenceFilteredFar.tsv"
    positions = set()

    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if int(row["ptm_pos"]) != int(ptm_pos):
                continue
            if uniprot and row.get("UniProt") != uniprot:
                continue

            fields = [row.get("near_mutation", ""), row.get("far_mutations_prevalence_filtered", "")]
            for field in fields:
                for token in field.split(","):
                    match = MUT_RE.search(token.strip())
                    if match:
                        positions.add(int(match.group(1)))

    return sorted(positions)


target_uniprot = "P35222" #change to be UniProt ID of interest (must match UniProt column in TSV)
ptm_positions = parse_ptm_positions(target_uniprot)

if not ptm_positions:
    raise ValueError(f"No PTM positions found in TSV for UniProt {target_uniprot}.")

for ptm_position in ptm_positions:
    mutation_positions = parse_mutation_positions(ptm_position, uniprot=target_uniprot)
    nearby = find_nearby_mutations(chain, ptm_position, mutation_positions)

    print(f"\nPTM position {ptm_position} ({target_uniprot}):")
    if nearby:
        for hit in nearby:
            print(f"  Mutation at {hit['mutation_pos']} is {hit['distance']:.2f} Å away")
    else:
        print("  No nearby mutations found within cutoff distance.")

from __future__ import annotations

import argparse
import os
import re
import time
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

import pandas as pd
import requests

AF_PRED_ENDPOINT = "https://alphafold.ebi.ac.uk/api/prediction/{}"
ACC_RE = re.compile(r"^(?:[A-NR-Z0-9][A-Z0-9]{5}|[OPQ][0-9][A-Z0-9]{3}[0-9])(?:-\d+)?$")


def clean_accession(cell: Any) -> Optional[str]:
    if cell is None:
        return None
    s = str(cell).strip()
    if not s:
        return None
    parts = re.split(r"[;,\s|]+", s)
    for p in parts:
        p = p.strip()
        if ACC_RE.match(p):
            return p
    return None


def pick_urls(record: dict, prefer: str = "cif") -> dict[str, str]:
    urls = [v for v in record.values() if isinstance(v, str) and v.startswith("http")]

    def pick_by_priority(priorities):
        best_url = ""
        best_score = -1
        for u in urls:
            lu = u.lower()
            for endings, score in priorities:
                if any(lu.endswith(e) for e in endings):
                    s = score + (1 if lu.endswith(".gz") else 0)
                    if s > best_score:
                        best_score = s
                        best_url = u
        return best_url

    if prefer == "cif":
        structure_url = pick_by_priority([
            ((".cif.gz", ".cif"), 300),
            ((".bcif.gz", ".bcif"), 200),
            ((".pdb.gz", ".pdb"), 100),
        ])
    else:
        structure_url = pick_by_priority([
            ((".pdb.gz", ".pdb"), 300),
            ((".cif.gz", ".cif"), 200),
            ((".bcif.gz", ".bcif"), 100),
        ])
    pae_url = pick_by_priority([
        ((".pae.json.gz", ".pae.json"), 300),
        ((".json.gz", ".json"), 100),
    ])
    return {"structure_url": structure_url, "pae_url": pae_url}

def fetch_prediction(acc: str, session: requests.Session, retries: int = 4) -> Optional[Any]:
    url = AF_PRED_ENDPOINT.format(acc)
    backoff = 1.6
    for attempt in range(retries):
        r = session.get(url, timeout=30)
        if r.status_code == 200:
            return r.json()
        if r.status_code == 404:
            return None
        time.sleep(backoff ** attempt)
    r.raise_for_status()
    return None

def download(url: str, outpath: Path, session: requests.Session, retries: int = 4) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    if outpath.exists() and outpath.stat().st_size > 0:
        return
    backoff = 1.6
    for attempt in range(retries):
        with session.get(url, stream=True, timeout=90) as r:
            if r.status_code == 200:
                tmp = outpath.with_suffix(outpath.suffix + ".part")
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 256):
                        if chunk:
                            f.write(chunk)
                os.replace(tmp, outpath)
                return
        time.sleep(backoff ** attempt)
    raise RuntimeError(f"Failed download after retries: {url}")


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p)
    # default TSV/CSV guess
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    return pd.read_csv(p, sep="\t")


def main(in_path: str, id_column: str, out_dir: str, prefer: str, also_pae: bool, delay: float) -> None:
    df = read_table(in_path)
    if id_column not in df.columns:
        raise ValueError(f"Column '{id_column}' not found. Columns: {list(df.columns)}")

    raw = df[id_column].tolist()
    accs = sorted({a for a in (clean_accession(x) for x in raw) if a})

    out_base = Path(out_dir)
    out_base.mkdir(parents=True, exist_ok=True)

    report_rows = []

    with requests.Session() as s:
        s.headers.update({"User-Agent": "bulk-afdb-downloader/1.0"})
        for i, acc in enumerate(accs, 1):
            row = {"UniProt": acc, "status": "", "structure_file": "", "pae_file": "", "note": ""}
            try:
                meta = fetch_prediction(acc, s)
                if not meta:
                    row["status"] = "NO_ENTRY"
                    row["note"] = "No AlphaFold DB record (404)"
                    report_rows.append(row)
                    print(f"[{i}/{len(accs)}] {acc}: no entry")
                    continue

                record = meta[0] if isinstance(meta, list) and meta else meta
                urls = pick_urls(record, prefer=prefer)

                if not urls.get("structure_url"):
                    row["status"] = "NO_STRUCTURE_URL"
                    row["note"] = f"No structure URL found; keys={list(record.keys())[:12]}..."
                    report_rows.append(row)
                    print(f"[{i}/{len(accs)}] {acc}: no structure url")
                    continue

                struct_url = urls["structure_url"]
                struct_name = struct_url.split("/")[-1]
                struct_path = out_base / acc / struct_name
                download(struct_url, struct_path, s)
                row["structure_file"] = str(struct_path)
                row["status"] = "DOWNLOADED"

                if also_pae and urls.get("pae_url"):
                    pae_url = urls["pae_url"]
                    if "pae" in pae_url.lower():
                        pae_name = pae_url.split("/")[-1]
                        pae_path = out_base / acc / pae_name
                        download(pae_url, pae_path, s)
                        row["pae_file"] = str(pae_path)

                report_rows.append(row)
                print(f"[{i}/{len(accs)}] {acc}: ok -> {struct_path.name}")
                time.sleep(delay)

            except Exception as e:
                row["status"] = "ERROR"
                row["note"] = str(e)
                report_rows.append(row)
                print(f"[{i}/{len(accs)}] {acc}: ERROR {e}")

    rep = pd.DataFrame(report_rows)
    rep_path = out_base / "download_report.tsv"
    rep.to_csv(rep_path, sep="\t", index=False)
    print(f"\nWrote report: {rep_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("in_path", help="Input .tsv/.csv/.xlsx file")
    ap.add_argument("--id_column", default="UniProt", help="Column containing UniProt accessions")
    ap.add_argument("--out_dir", default="afdb_models", help="Output directory")
    ap.add_argument("--prefer", choices=["cif", "pdb"], default="cif", help="Prefer mmCIF or PDB when both available")
    ap.add_argument("--also_pae", action="store_true", help="Also download PAE json when url contains 'pae'")
    ap.add_argument("--delay", type=float, default=0.1, help="Polite delay between requests (seconds)")
    args = ap.parse_args()

    main(args.in_path, args.id_column, args.out_dir, args.prefer, args.also_pae, args.delay)

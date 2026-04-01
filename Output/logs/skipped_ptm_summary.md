# Skipped PTM Summary

49 proteins had PTMs excluded from analysis for non-biological reasons.
Each entry lists the skip reason, PTM count, and what to do if you want to recover this data.

---

## No AFDB Entry (16 proteins)
These proteins have no AlphaFold DB record at all. All are very large proteins likely excluded during AFDB database construction. To recover: generate structures via the AlphaFold Server (https://alphafoldserver.com, limit 5000 AA) or ESMFold, place CIF in `cif_models/{UniProt}/AF-{UniProt}-F1-model_v6.cif`.

| UniProt | Gene | Skipped PTMs | Notes |
|---------|------|-------------|-------|
| Q8WXI7 | MUC16 | 278 | 14507 AA — too large for AlphaFold Server |
| P25054 | APC | 186 | 2843 AA — major tumor suppressor, worth modeling |
| Q96T58 | SPEN | 119 | 3664 AA |
| Q03164 | KMT2A | 110 | 3969 AA — recurrent leukemia fusion gene |
| Q13315 | ATM | 132 | 3056 AA — major DNA damage checkpoint kinase, worth modeling |
| P51587 | BRCA2 | 183 | 3418 AA — major breast/ovarian cancer suppressor, worth modeling |
| O14686 | KMT2D | 127 | 5537 AA — too large for AlphaFold Server |
| Q9Y4A5 | TRRAP | 60 | 3859 AA |
| Q8NEZ4 | KMT2C | 76 | 4911 AA — too large for AlphaFold Server |
| Q9NR09 | BIRC6 | 54 | 4857 AA — too large for AlphaFold Server |
| O95071 | UBR5 | 85 | 2799 AA |
| Q14517 | FAT1 | 76 | 4588 AA — too large for AlphaFold Server |
| Q15911 | ZFHX3 | 38 | 3703 AA |
| Q6V0I7 | FAT4 | 13 | 4981 AA — too large for AlphaFold Server |
| P49792 | RANBP2 | 306 | 3224 AA |
| Q9NZR2 | LRP1B | 8 | 4599 AA — too large for AlphaFold Server |

---

## No Canonical AFDB Model (7 proteins)
AFDB has entries for these proteins but only modeled specific isoforms — never the canonical sequence. All are large proteins. To recover: same approach as above, model the canonical sequence manually.

| UniProt | Gene | Skipped PTMs | Notes |
|---------|------|-------------|-------|
| Q63HN8 | RNF213 | 100 | 5207 AA canonical — too large for AlphaFold Server; AFDB only has 557 AA and 1063 AA isoforms |
| P21359 | NF1 | 73 | 2839 AA — major RAS pathway tumor suppressor, worth modeling |
| Q8NFP9 | NBEA | 20 | 2946 AA — AFDB has 2 isoform models only |
| Q7Z407 | CSMD3 | 4 | 3707 AA |
| Q99102 | MUC4 | 4 | 5412 AA — too large for AlphaFold Server |
| Q99996 | AKAP9 | 2 | 3907 AA |
| Q8TDW7 | FAT3 | 1 | 4557 AA — too large for AlphaFold Server |

---

## Residue Mismatch (25 proteins)
The amino acid at the PTM position in the canonical AlphaFold structure does not match the residue in the PTMD annotation. Causes: PTMD annotated against an outdated UniProt sequence version, an alternatively spliced exon in the source isoform, or the entire annotation being from a different isoform. PTMD does not record which sequence version or isoform it used. These PTMs cannot be mapped to the canonical structure without resolving the source sequence.

| UniProt | Gene | Skipped PTMs | Likely cause |
|---------|------|-------------|--------------|
| Q8NFD5 | ARID1B | 33 | 34/35 PTMs mismatch — entire annotation on wrong isoform |
| P52948 | NUP98 | 20 | 20/26 mismatch — likely isoform-specific region |
| P78549 | NTHL1 | 10 | 10/10 mismatch — entire annotation on wrong isoform |
| Q96RK0 | CIC | 11 | 11/68 mismatch — scattered, likely isoform-specific exons |
| Q9HB09 | BCL2L12 | 9 | 9/9 mismatch — entire annotation on wrong isoform |
| O00255 | MEN1 | 5 | mixed near-end and internal — likely partially isoform-specific |
| O75030 | MITF | 5 | scattered mismatches — isoform-specific |
| Q6NWY9 | PRPF40B | 5 | 5/6 mismatch, clustered near C-terminus — alternative C-terminal exon |
| E9PAV3 | NACA | 4 | clustered at N-terminus (positions 43–191) — alternative N-terminal exon |
| P55197 | MLLT10 | 4 | 4 consecutive mismatches (696–709) — alternatively spliced exon |
| P19544 | WT1 | 3 | also has position_not_in_structure; isoform-specific |
| Q9HBE5 | IL21R | 3 | 3/3 mismatch — entire annotation on wrong isoform |
| P38936 | CDKN1A | 7 | also has position_not_in_structure; likely wrong isoform throughout |
| P30622 | CLIP1 | 2 | 2 mismatches spread out — sequence version drift |
| P16220 | CREB1 | 2 | 2/2 mismatch — entire annotation on wrong isoform |
| P63092 | GNAS | 2 | also has position_not_in_structure; canonical is 394 AA, PTMD uses longer isoform |
| O15013 | ARHGEF10 | 3 | 2 mismatches near C-terminus — alternative C-terminal exon |
| P01106 | MYC | 1 | T58→P58 — likely sequence version update |
| P08575 | PTPRC | 1 | S973→R973 — likely sequence version update |
| P10275 | AR | 1 | K630→R630 — conservative substitution, sequence version update |
| P39880 | CUX1 | 1 | S613→F613 — sequence version update |
| Q03112 | MECOM | 1 | S538→T538 — conservative substitution, sequence version update |
| Q15910 | EZH2 | 1 | K27→R27 — conservative substitution, sequence version update |
| Q5H9F3 | BCORL1 | 1 | near C-terminus — alternative C-terminal exon |
| Q86Y26 | NUTM1 | 1 | near C-terminus — alternative C-terminal exon |

---

## Position Beyond Canonical Length (4 proteins)
PTM positions in PTMD exceed the canonical protein length entirely. These positions do not exist in any known sequence of the protein. Likely PTMD annotation errors or annotations from an undocumented longer isoform.

| UniProt | Gene | Skipped PTMs | Detail |
|---------|------|-------------|--------|
| P38936 | CDKN1A | 10 | Canonical length 164 AA; PTMs annotated at positions 175–197 — positions do not exist |
| P63092 | GNAS | 28 | Canonical length 394 AA; PTMs annotated at positions 532–694 — likely long isoform (XLas, ~820 AA) |
| P19544 | WT1 | 1 | Canonical length 449 AA; S461 annotated beyond end |
| P31749 | AKT1 | 1 | Canonical length 480 AA; T554 annotated beyond end |

---

## Priority Candidates for Manual Recovery

If the group wants to manually obtain structures and recover data, these are the highest-value targets:

1. **ATM (Q13315)** — 132 PTMs, 3056 AA, fits AlphaFold Server, major cancer gene
2. **BRCA2 (P51587)** — 183 PTMs, 3418 AA, fits AlphaFold Server, major cancer gene
3. **APC (P25054)** — 186 PTMs, 2843 AA, fits AlphaFold Server, major cancer gene
4. **NF1 (P21359)** — 73 PTMs, 2839 AA, fits AlphaFold Server, isoform-only in AFDB
5. **KMT2A (Q03164)** — 110 PTMs, 3969 AA, fits AlphaFold Server, recurrent leukemia gene

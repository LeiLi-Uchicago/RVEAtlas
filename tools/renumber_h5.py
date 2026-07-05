#!/usr/bin/env python3
"""Renumber an H5 HA biological-assembly-1 trimer PDB to the app's H5 numbering.

Why: unlike 4JTV (H1) and 4HMG (H3), whose PDB auth residue numbers already equal
the app's mature HA1/HA2 numbering, crystal H5 structures use H3-style numbering
with insertion codes (e.g. 2FK0 / 9NRR chain A: 125A/125B/133A/...). No constant
offset reconciles that with the app's sequential H5 numbering, so we align the app
H5 sequence to the structure and rewrite the auth residue numbers. After this,
structure_config.tsv can keep resi_from_numbering = TRUE.

Requirements of the input assembly: a single MODEL with distinct protein chains
A-F (HA1 = A,C,E; HA2 = B,D,F) plus optional glycan chains. Assemblies deposited
as several MODELs with repeated chain IDs are NOT supported (they'd collide).

Inputs (run from repo root):
  - <input.pdb1>                      RCSB biological assembly 1, e.g.
      curl -sL https://files.rcsb.org/download/2FK0.pdb1.gz | gunzip > 2FK0.pdb1
  - ha_numbering_review_table.csv     app numbering (H5NX Manual_Aligned_AA per Numbering_Position)

Output:
  - <output.pdb>  waters/buffer stripped, glycans kept; chains A/C/E (HA1) and
      B/D/F (HA2) renumbered so auth resSeq == app Numbering_Position (iCode cleared);
      structure-only insertions get sentinel numbers >=1000 so they never match an
      epitope/selected site. SEQRES/HELIX/SHEET/SSBOND/LINK dropped (3Dmol recomputes SS).

Usage:
  python3 tools/renumber_h5.py 2FK0.pdb1 www/structures/2FK0_trimer.pdb
  python3 tools/renumber_h5.py 9NRR.pdb1 www/structures/9NRR_trimer.pdb
"""
import csv
import sys

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
NUMBERING = 'ha_numbering_review_table.csv'
# HA1 -> chain A, HA2 -> chain B (protomer 1); the same numbering applies to C/E and D/F.
REGION_CHAIN = {'HA1': 'A', 'HA2': 'B'}
CHAIN_REGION = {'A': 'HA1', 'C': 'HA1', 'E': 'HA1', 'B': 'HA2', 'D': 'HA2', 'F': 'HA2'}
# Solvent / buffer / cryoprotectant / ion HETATMs to drop; glycans (NAG/BMA/MAN/FUC/
# GAL/SIA/NDG/...) and everything else are kept.
SOLVENT = {'HOH', 'EDO', 'PEG', 'PG4', 'PG0', 'GOL', 'SO4', 'PO4', 'ACT', 'DMS', 'MPD',
           'MES', 'EPE', 'TRS', 'FMT', 'NO3', 'IOD', 'BR', 'CL', 'NA', 'K', 'MG', 'CA',
           'ZN', 'MN', 'NI', 'CD', 'FLC', 'CIT', 'BME', 'IPA', 'PGE'}


def app_seq(region):
    rows = []
    for r in csv.DictReader(open(NUMBERING)):
        if r['Subtype'] != 'H5NX' or r['HA_Region'] != region:
            continue
        if str(r['Is_Alignment_Gap']).lower() in ('true', '1'):
            continue
        aa = r['Manual_Aligned_AA']
        if aa in AA3.values():
            rows.append((int(float(r['Numbering_Position'])), aa))
    rows.sort()
    return rows


def struct_res(pdb, chain):
    out, inmodel = [], 0
    for line in open(pdb):
        if line.startswith('MODEL'):
            inmodel += 1
        if inmodel > 1:
            break
        if line.startswith('ATOM') and line[21] == chain and line[12:16] == ' CA ':
            out.append((int(line[22:26]), line[26], AA3.get(line[17:20].strip(), 'X')))
    return out


def needleman_wunsch(a, b, match=2, mis=-1, gap=-2):
    n, m = len(a), len(b)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if a[i - 1][1] == b[j - 1][2] else mis
            F[i][j] = max(F[i - 1][j - 1] + s, F[i - 1][j] + gap, F[i][j - 1] + gap)
    i, j, ali = n, m, []
    while i > 0 and j > 0:
        s = match if a[i - 1][1] == b[j - 1][2] else mis
        if F[i][j] == F[i - 1][j - 1] + s:
            ali.append((i - 1, j - 1)); i -= 1; j -= 1
        elif F[i][j] == F[i - 1][j] + gap:
            ali.append((i - 1, None)); i -= 1
        else:
            ali.append((None, j - 1)); j -= 1
    while i > 0:
        ali.append((i - 1, None)); i -= 1
    while j > 0:
        ali.append((None, j - 1)); j -= 1
    ali.reverse()
    return ali


def build_maps(pdb):
    maps = {}
    for region, chain in REGION_CHAIN.items():
        a, b = app_seq(region), struct_res(pdb, chain)
        ali = needleman_wunsch(a, b)
        mp, sentinel, matched, mism = {}, 1000, 0, 0
        for ai, bj in ali:
            if bj is None:
                continue
            key = (b[bj][0], b[bj][1])
            if ai is None:
                mp[key] = sentinel; sentinel += 1
            else:
                mp[key] = a[ai][0]; matched += 1; mism += (a[ai][1] != b[bj][2])
        maps[region] = mp
        print(f"{region}: matched {matched}/{len(b)} struct residues; {mism} AA mismatches; app residues={len(a)}")
    return maps


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: python3 tools/renumber_h5.py <input.pdb1> <output.pdb>")
    pdb, out_path = sys.argv[1], sys.argv[2]
    maps = build_maps(pdb)
    inmodel = 0
    with open(out_path, 'w') as out:
        for line in open(pdb):
            rec = line[:6]
            if line.startswith('MODEL'):
                inmodel += 1
                if inmodel > 1:
                    break
                continue
            if line.startswith('ENDMDL'):
                continue
            if rec == 'HETATM' and line[17:20].strip() in SOLVENT:
                continue
            if rec in ('ATOM  ', 'HETATM'):
                ch = line[21]
                if rec == 'ATOM  ' and ch in CHAIN_REGION:
                    key = (int(line[22:26]), line[26])
                    nn = maps[CHAIN_REGION[ch]].get(key, 9000 + int(line[22:26]) % 900)
                    line = line[:22] + f"{nn:4d}" + " " + line[27:]  # resSeq cols 23-26, clear iCode col 27
                out.write(line)
            elif rec == 'TER   ':
                out.write(line)
    print(f"wrote {out_path}")


if __name__ == '__main__':
    main()

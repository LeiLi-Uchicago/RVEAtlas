#!/usr/bin/env python3
"""Renumber an influenza NA biological-assembly tetramer PDB to the app's NA numbering.

NA is a homotetramer whose RCSB biological assembly is usually delivered as several
MODELs (or one MODEL with a single chain plus 4-fold symmetry) rather than four
distinct chains in one model. This script flattens the assembly: every protein
protomer becomes a distinct chain (A, B, C, D, ...) in a single MODEL, and each is
renumbered so its auth resSeq equals the app's sequential NA position. The app
position is obtained by globally aligning the app NA consensus (a FASTA) to the
protomer sequence, so the structure's own author numbering is irrelevant.

After this, structure_config.tsv can use resi_from_numbering = FALSE (identity:
app position == PDB resi) with region_chains "NA:A,B,C,D".

Inputs (run from repo root):
  - <input.pdb1>   RCSB biological assembly 1, e.g.
      curl -sL https://files.rcsb.org/download/3NSS.pdb1.gz | gunzip > 3NSS.pdb1
  - <app_seq.fa>   FASTA whose first record is the app NA consensus for this subtype

Output:
  - <output.pdb>   waters/buffer stripped; glycans and the catalytic Ca2+ kept;
      protomers renamed A,B,C,D and renumbered to app positions (iCode cleared);
      structure-only residues that don't align get sentinel numbers >=100000 so they
      never collide with an epitope/selected site. SEQRES/HELIX/SHEET/SSBOND/LINK
      dropped (3Dmol recomputes secondary structure).

Usage:
  python3 tools/renumber_na.py 3NSS.pdb1 na_n1_app.fa www/structures/3NSS_tetramer.pdb
  python3 tools/renumber_na.py 4GZO.pdb1 na_n2_app.fa www/structures/4GZO_tetramer.pdb
"""
import sys

try:
    from Bio import pairwise2
except ImportError:
    sys.exit("biopython required: pip install biopython")

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
# Solvent / buffer / cryoprotectant to drop. Glycans (NAG/BMA/MAN/FUC/GAL/SIA/...)
# and the catalytic calcium (CA) are kept.
SOLVENT = {'HOH', 'EDO', 'PEG', 'PG4', 'PG0', 'GOL', 'SO4', 'PO4', 'ACT', 'DMS', 'MPD',
           'MES', 'EPE', 'TRS', 'FMT', 'NO3', 'IOD', 'BR', 'CL', 'NA', 'K', 'MG',
           'ZN', 'MN', 'NI', 'CD', 'FLC', 'CIT', 'BME', 'IPA', 'PGE', 'HEPES'}
NEW_CHAINS = "ABCDEFGHIJKLMNOP"


def read_app_seq(path):
    seq = []
    started = False
    for line in open(path):
        if line.startswith('>'):
            if started:
                break
            started = True
            continue
        seq.append(line.strip())
    return "".join(seq)


def protomers(lines):
    """Yield (proto_index, list_of_record_lines) for each protein protomer, in
    order of appearance across MODELs. A protomer = one (model, chain) that has
    ATOM protein records."""
    model = 1
    buckets = {}
    order = []
    for line in lines:
        rec = line[:6]
        if rec == 'MODEL ':
            model = int(line[10:14])
            continue
        if rec in ('ATOM  ', 'HETATM'):
            chain = line[21]
            key = (model, chain)
            if key not in buckets:
                buckets[key] = []
                order.append(key)
            buckets[key].append(line)
    for key in order:
        recs = buckets[key]
        has_protein = any(r[:6] == 'ATOM  ' and r[17:20] in AA3 for r in recs)
        if has_protein:
            yield key, recs


def residue_seq(recs):
    """Ordered list of (resSeq, iCode, aa1) for the protein residues of a protomer."""
    out = []
    seen = set()
    for r in recs:
        if r[:6] != 'ATOM  ':
            continue
        resname = r[17:20]
        if resname not in AA3:
            continue
        rid = (r[22:26], r[26])
        if rid in seen:
            continue
        seen.add(rid)
        out.append((r[22:26], r[26], AA3[resname]))
    return out


def build_map(app_seq, res_list):
    """Align app_seq to the protomer sequence; return dict (resSeq,iCode)->app_pos."""
    struct_seq = "".join(aa for _, _, aa in res_list)
    aln = pairwise2.align.globalms(app_seq, struct_seq, 2, -1, -10, -0.5,
                                   one_alignment_only=True)[0]
    a_app, a_str = aln.seqA, aln.seqB
    i_app = i_str = 0
    mapping = {}
    for c_app, c_str in zip(a_app, a_str):
        if c_app != '-':
            i_app += 1
        if c_str != '-':
            if c_str != '-' and c_app != '-':
                resSeq, icode, _ = res_list[i_str]
                mapping[(resSeq, icode)] = i_app
            i_str += 1
    return mapping


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    inp, appfa, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    app_seq = read_app_seq(appfa)
    lines = open(inp).read().splitlines()

    out = []
    # Sentinel base for structure-only (unalignable) residues. Kept above any real
    # protein length so sv_structure_resolved_positions() can drop them without
    # hiding genuine high-numbered residues (e.g. Spike goes to 1273).
    sentinel = 100000
    matched_total = struct_total = 0
    for idx, ((model, chain), recs) in enumerate(protomers(lines)):
        if idx >= len(NEW_CHAINS):
            break
        new_chain = NEW_CHAINS[idx]
        res_list = residue_seq(recs)
        mapping = build_map(app_seq, res_list)
        struct_total += len(res_list)
        matched_total += len(mapping)
        for r in recs:
            rec = r[:6]
            resname = r[17:20]
            if rec == 'HETATM':
                if resname.strip() in SOLVENT:
                    continue  # drop waters/buffer
                # keep glycans + Ca2+; retain their own numbering, remap chain
                out.append(r[:21] + new_chain + r[22:])
                continue
            if rec != 'ATOM  ':
                continue
            if resname not in AA3:
                continue
            rid = (r[22:26], r[26])
            app_pos = mapping.get(rid)
            if app_pos is None:
                app_pos = sentinel
                sentinel += 1
            out.append(r[:21] + new_chain + ("%4d" % app_pos) + ' ' + r[27:])

    with open(outp, 'w') as fh:
        fh.write("REMARK  Renumbered to app NA numbering by tools/renumber_na.py\n")
        for line in out:
            fh.write(line.rstrip('\n') + "\n")
        fh.write("END\n")

    print("protomers: %d; residues matched %d/%d"
          % (min(len(NEW_CHAINS), sum(1 for _ in protomers(lines))),
             matched_total, struct_total))


if __name__ == '__main__':
    main()

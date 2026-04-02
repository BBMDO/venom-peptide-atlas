#!/usr/bin/env python3
"""
param_alphafold_updated.py
Structure-based analytics for PDB/mmCIF (incl. AlphaFold models).
"""
import os, argparse, numpy as np, pandas as pd
from glob import glob
from typing import Dict, Tuple, List, Optional
from Bio.PDB import PDBParser, MMCIFParser, DSSP, Polypeptide
from Bio.PDB.Polypeptide import is_aa
import freesasa

EISENBERG_SCALE = {'A':0.62,'R':-2.53,'N':-0.78,'D':-0.90,'C':0.29,'Q':-0.85,'E':-0.74,'G':0.48,'H':-0.40,'I':1.38,'L':1.06,'K':-1.50,'M':0.64,'F':1.19,'P':0.12,'S':-0.18,'T':-0.05,'W':0.81,'Y':0.26,'V':1.08}
THREE2ONE = {'ALA':'A','VAL':'V','ILE':'I','LEU':'L','MET':'M','PHE':'F','TRP':'W','TYR':'Y','CYS':'C','GLY':'G','PRO':'P','SER':'S','THR':'T','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R'}

def get_parser(path: str):
    ext = os.path.splitext(path)[1].lower()
    return MMCIFParser(QUIET=True) if ext in (".cif",".mmcif") else PDBParser(QUIET=True)
def load_structure(path: str):
    return get_parser(path).get_structure("X", path)
def chain_sequence(chain) -> str:
    ppb = Polypeptide.PPBuilder()
    seqs = [str(pp.get_sequence()) for pp in ppb.build_peptides(chain)]
    return "".join(seqs)
def dssp_secfrac(structure, model_id=0):
    try:
        dssp = DSSP(structure[model_id], None)
        ss = [s[2] for s in dssp]; n = len(ss) or 1
        helix = sum(x in ('H','G','I') for x in ss)/n
        sheet = sum(x in ('E','B') for x in ss)/n
        coil = 1.0 - helix - sheet
        return {"sec_helix":helix,"sec_sheet":sheet,"sec_coil":coil, "dssp_len": n}
    except Exception:
        return {"sec_helix":np.nan,"sec_sheet":np.nan,"sec_coil":np.nan, "dssp_len": np.nan}
def freesasa_features(path: str):
    try:
        s = freesasa.Structure(path); r = freesasa.calc(s); sasa_total = r.totalArea()
        hydrophobic = set(list("AVILMFWYC"))
        res_areas = {}
        for i in range(s.nAtoms()):
            resn = s.residueName(i).strip().upper(); chain = s.chainLabel(i); resi = s.residueNumber(i)
            key = (chain,resi,resn); res_areas[key] = res_areas.get(key,0.0) + r.atomArea(i)
        hyd_sasa = 0.0
        for (chain,resi,resn), area in res_areas.items():
            one = THREE2ONE.get(resn, None)
            if one in hydrophobic: hyd_sasa += area
        frac_hyd_surf = hyd_sasa/sasa_total if sasa_total>0 else np.nan
        return {"sasa_total":sasa_total, "frac_hydrophobic_surface": frac_hyd_surf}
    except Exception:
        return {"sasa_total":np.nan, "frac_hydrophobic_surface": np.nan}
def radius_of_gyration(structure, model_id=0, chain_id=None):
    try:
        coords = []
        model = structure[model_id]
        for ch in model:
            if chain_id and ch.id != chain_id: continue
            for res in ch:
                if is_aa(res, standard=True) and 'CA' in res:
                    coords.append(res['CA'].get_coord())
        if not coords: return np.nan
        X = np.vstack(coords); center = X.mean(axis=0)
        rg2 = ((X - center)**2).sum(axis=1).mean()
        return float(np.sqrt(rg2))
    except Exception:
        return np.nan
def count_disulfide(structure, model_id=0, cutoff=2.2):
    try:
        model = structure[model_id]; sg = []
        for ch in model:
            for res in ch:
                if is_aa(res, standard=True) and res.get_resname()=="CYS" and 'SG' in res:
                    sg.append(res['SG'].get_coord())
        cnt = 0
        for i in range(len(sg)):
            for j in range(i+1,len(sg)):
                if np.linalg.norm(sg[i]-sg[j]) <= cutoff: cnt += 1
        return cnt
    except Exception:
        return 0
def extract_plddt(chain):
    try:
        vals = [atom.get_bfactor() for res in chain for atom in res]
        if not vals: return (np.nan,np.nan,np.nan)
        arr = np.array(vals, dtype=float)
        return (float(np.mean(arr)), float(np.median(arr)), float(np.mean(arr<70.0)))
    except Exception:
        return (np.nan,np.nan,np.nan)
def contact_density(structure, model_id=0, cutoff=8.0):
    try:
        coords = []
        model = structure[model_id]
        for ch in model:
            for res in ch:
                if is_aa(res, standard=True) and 'CA' in res:
                    coords.append(res['CA'].get_coord())
        n = len(coords)
        if n < 2: return np.nan
        X = np.vstack(coords); cnt = 0
        for i in range(n):
            for j in range(i+1,n):
                if np.linalg.norm(X[i]-X[j]) < cutoff: cnt += 1
        return float(2.0*cnt/n)
    except Exception:
        return np.nan
def process_path(path: str) -> pd.DataFrame:
    if os.path.isdir(path):
        files = sorted(glob(os.path.join(path, "*.pdb")) + glob(os.path.join(path, "*.cif")) + glob(os.path.join(path, "*.mmcif")))
    else:
        files = [path]
    rows = []
    for fp in files:
        try:
            structure = load_structure(fp)
        except Exception:
            rows.append({"input_file": os.path.basename(fp), "error":"parse_failed"}); continue
        model = next(structure.get_models())
        for chain in model:
            if not any(is_aa(r, standard=True) for r in chain): continue
            seq = chain_sequence(chain)
            sec = dssp_secfrac(structure, model.id)
            sasa = freesasa_features(fp)
            rg = radius_of_gyration(structure, model.id, chain.id)
            ssc = count_disulfide(structure, model.id)
            plddt_mean, plddt_median, frac_plddt_lt70 = extract_plddt(chain)
            cden = contact_density(structure, model.id)
            rows.append({"input_file": os.path.basename(fp), "chain": chain.id, "n_residues": len(seq), **sec, **sasa, "rg": rg, "disulfide_pairs": ssc, "plddt_mean": plddt_mean, "plddt_median": plddt_median, "frac_plddt_lt70": frac_plddt_lt70, "contact_density": cden})
    return pd.DataFrame(rows)
def main():
    ap = argparse.ArgumentParser(description="Structure-based peptide analytics (PDB/mmCIF).")
    ap.add_argument("-p","--path", required=True, help="PDB/mmCIF file or directory")
    ap.add_argument("-o","--output", required=True, help="Output CSV")
    args = ap.parse_args()
    df = process_path(args.path); df.to_csv(args.output, index=False)
    print(f"[OK] Wrote {args.output} with {len(df)} rows.")
if __name__ == "__main__":
    main()

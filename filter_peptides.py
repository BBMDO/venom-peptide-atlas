#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filtra e limpa peptídeos de um FASTA, calcula propriedades
e gera:
  - FASTA com sequências aprovadas (apenas limpeza/validação)
  - TSV 'aprovados' com propriedades
  - TSV 'descartes' apenas para itens inválidos (vazios ou erro de cálculo)

IMPORTANTE:
- Não há descarte por limiares de propriedades (pI/GRAVY/carga/instabilidade).
- Apenas sequências que ficam vazias após a limpeza são descartadas,
  ou registros que gerem erro durante o cálculo de propriedades.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Alfabeto de 20 aa canônicos; mapeamos U (Sec) -> C e removemos o resto
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
PH_FOR_CHARGE = 7.0  # pH usado para calcular carga

def clean_sequence(seq_str: str) -> str:
    """Normaliza para maiúsculas, remove *, espaços/quebras e caracteres não canônicos.
       Mapeia U (selenocisteína) -> C (ajuste se não desejar)."""
    s = str(seq_str).upper()
    s = s.replace("*", "").replace(" ", "").replace("\r", "").replace("\n", "")
    s = s.replace("U", "C")  # opcional
    s = "".join(ch for ch in s if ch in VALID_AA)
    return s

def calculate_properties(seq_clean: str):
    """Calcula (pI, GRAVY, charge@pH, instability) usando Biopython."""
    analysed = ProteinAnalysis(seq_clean)
    pI = analysed.isoelectric_point()
    gravy = analysed.gravy()
    charge = analysed.charge_at_pH(PH_FOR_CHARGE)
    instability = analysed.instability_index()
    return pI, gravy, charge, instability

def main():
    parser = argparse.ArgumentParser(
        description="Limpa sequências, calcula propriedades e salva aprovados (sem aplicar limiares)."
    )
    parser.add_argument("-i", "--input", required=True, help="FASTA de entrada")
    parser.add_argument("-o", "--output", required=True, help="FASTA de saída com sequências aprovadas")
    parser.add_argument("-l", "--log", required=True, help="TSV de descartes (motivos)")
    parser.add_argument("-r", "--report", required=True, help="TSV com propriedades dos aprovados")
    args = parser.parse_args()

    input_fa = args.input
    output_fasta = args.output
    log_tsv = args.log
    report_tsv = args.report

    total = 0
    approved_records = []
    # Tabelas
    discarded_rows = [["id", "reason", "length_raw", "length_clean", "note"]]
    approved_rows = [["id", "length_clean", "pI", "GRAVY", f"charge_pH{PH_FOR_CHARGE}", "instability"]]

    for record in SeqIO.parse(input_fa, "fasta"):
        total += 1
        raw_seq = str(record.seq)
        seq_clean = clean_sequence(raw_seq)

        # Descarte 1: vazia após limpeza
        if len(seq_clean) == 0:
            discarded_rows.append([record.id, "empty_after_cleaning", len(raw_seq), 0,
                                   "sequência ficou vazia após remoção de caracteres não canônicos/*/espaços"])
            continue

        # Calcula propriedades (sem aplicar limiares de exclusão)
        try:
            pI, gravy, charge, instability = calculate_properties(seq_clean)
        except Exception as e:
            discarded_rows.append([record.id, "calc_error", len(raw_seq), len(seq_clean),
                                   f"erro no cálculo de propriedades: {e}"])
            continue

        # Aprovado (sempre que passou pela limpeza e cálculo)
        approved_records.append(SeqRecord(Seq(seq_clean), id=record.id, description=record.description))
        approved_rows.append([
            record.id, len(seq_clean),
            f"{pI:.4f}", f"{gravy:.4f}", f"{charge:.4f}", f"{instability:.4f}"
        ])

    # Grava FASTA dos aprovados
    with open(output_fasta, "w") as fh_fa:
        SeqIO.write(approved_records, fh_fa, "fasta")

    # Grava TSV dos aprovados
    with open(report_tsv, "w", encoding="utf-8") as fh_ok:
        for row in approved_rows:
            fh_ok.write("\t".join(map(str, row)) + "\n")

    # Grava TSV de descartes
    with open(log_tsv, "w", encoding="utf-8") as fh_log:
        for row in discarded_rows:
            fh_log.write("\t".join(map(str, row)) + "\n")

    # Resumo no STDERR
    print(f"[i] Total lidas: {total}", file=sys.stderr)
    print(f"[i] Aprovadas: {len(approved_records)}  → {output_fasta}", file=sys.stderr)
    print(f"[i] Propriedades (aprovados): {len(approved_rows)-1}  → {report_tsv}", file=sys.stderr)
    print(f"[i] Descartadas: {len(discarded_rows)-1}  → {log_tsv}", file=sys.stderr)

if __name__ == "__main__":
    main()

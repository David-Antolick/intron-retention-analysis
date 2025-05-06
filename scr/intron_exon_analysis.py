#!/usr/bin/env python3
import pandas as pd
from collections import defaultdict
from Bio.Seq import Seq
import math

def main():
    # Load Excel file and read sheets
    xls = pd.ExcelFile(input_excel)
    exon_df = pd.read_excel(xls, sheet_name=exon_sheet)
    exon_df.columns = exon_df.columns.str.strip()
    master_df = pd.read_excel(xls, sheet_name=master_sheet)
    master_df.columns = master_df.columns.str.strip()
    master_df["Ensembl ID"] = master_df["Ensembl ID"].astype(str).str.strip()

    # Build dictionary mapping Ensembl ID to its exon annotations
    exons_by_gene = {}
    for _, row in exon_df.iterrows():
        ensembl_id = row["ensembl_id"]
        if ensembl_id not in exons_by_gene:
            exons_by_gene[ensembl_id] = []
        exons_by_gene[ensembl_id].append({
            "start": row["coding_start"],
            "end": row["coding_end"],
        })

    updated_tps = {}
    for sheet in timepoint_sheets:
        if sheet in xls.sheet_names:
            tp_df = pd.read_excel(xls, sheet_name=sheet)
            tp_df.columns = tp_df.columns.str.strip()
            tp_df["Ensembl ID"] = tp_df["Ensembl ID"].astype(str).str.strip()
            
            # Merge with master to get gene sequence and coordinate data.
            needed_cols = ["Ensembl ID", "Gene Sequence", "gene_start", "gene_end", "strand"]
            tp_df = tp_df.merge(master_df[needed_cols], on="Ensembl ID", how="left", suffixes=("", "_master"))
            tp_df = tp_df.drop_duplicates(subset=["Ensembl ID", "intron_start", "intron_end"], keep="first")
            
            ir_full_seqs = []
            ir_proteins = []
            stop_to_intron_end_list = []
            stop_to_exon_junction_list = []
            
            for _, row in tp_df.iterrows():
                ensembl_id = row["Ensembl ID"]
                # If missing exon info or gene_start, skip processing.
                if ensembl_id not in exons_by_gene or pd.isna(row.get("gene_start")):
                    ir_full_seqs.append("")
                    ir_proteins.append("")
                    stop_to_intron_end_list.append(None)
                    stop_to_exon_junction_list.append(None)
                    continue
                
                # Reconstruct the intron-retained (IR) transcript.
                # This function returns the IR sequence, the absolute gene start (from the Desc field), 
                # and the nucleotide index (within the IR transcript) where the intron ends.
                ir_seq, gene_start, ir_intron_end = _reconstruct_intron_retained(row, exons_by_gene[ensembl_id])
                
                # Adjust for negative strand: reverse complement the IR transcript.
                strand = row.get("strand", "+")
                if strand in ("-", "-1", -1):
                    ir_seq = str(Seq(ir_seq).reverse_complement())
                    # NOTE: For negative strand genes, the IR transcript index for the intron end
                    # should be adjusted accordingly. For now, we use the same value.
                
                ir_full_seqs.append(ir_seq)
                
                # Translate the IR transcript starting at the first ATG.
                protein_seq, atg_index = _translate_dna_modified(ir_seq)
                ir_proteins.append(protein_seq)
                
                # Compute distances.
                # (a) Distance (in amino acids) from the first stop codon in the IR transcript to the end of the retained intron.
                # (b) Distance (in amino acids) from the first stop codon to the canonical next exon junction.
                dist_ir, canon_dist = find_stop_codon_distance(
                    protein_seq, 
                    atg_index if atg_index is not None else 0, 
                    ir_intron_end, 
                    gene_start, 
                    exons_by_gene[ensembl_id],
                    strand
                )
                stop_to_intron_end_list.append(dist_ir)
                stop_to_exon_junction_list.append(canon_dist)
            
            tp_df["ir_full_seq"] = ir_full_seqs
            tp_df["ir_protein"] = ir_proteins
            tp_df["stop_to_intron_end"] = stop_to_intron_end_list
            tp_df["stop_to_exon_junction"] = stop_to_exon_junction_list
            
            updated_tps[sheet] = tp_df
        else:
            updated_tps[sheet] = pd.DataFrame()

    # Write the updated data to a new Excel file.
    with pd.ExcelWriter(output_excel) as writer:
        exon_df.to_excel(writer, sheet_name=exon_sheet, index=False)
        master_df.to_excel(writer, sheet_name=master_sheet, index=False)
        for sheet in timepoint_sheets:
            df = updated_tps[sheet]
            if not df.empty:
                df.to_excel(writer, sheet_name=sheet, index=False)

def _reconstruct_intron_retained(row, exons):
    """
    Reconstruct the gene sequence with the intron retained.
    Assumes that:
      - The gene sequence (gs) is a full genomic sequence for the gene.
      - 'gene_start' is the absolute genomic coordinate of the gene’s start.
      - 'intron_start' and 'intron_end' (from the row) are provided as nucleotide offsets 
        relative to gene_start.
      - 'exons' is a list of exon dicts with absolute 'start' and 'end' coordinates.
    
    Returns:
      (ir_seq, gene_start, ir_intron_end)
      where ir_seq is the constructed IR transcript,
            gene_start is as provided,
            and ir_intron_end is the nucleotide index in the IR transcript where the retained intron ends.
    """
    gs = row.get("Gene Sequence", "")
    gs = "" if pd.isna(gs) else str(gs)
    intron_seq = row.get("Intron Sequence", "")
    if pd.isna(intron_seq):
        intron_seq = ""
    
    gene_start = row.get("gene_start", math.nan)
    intron_start = row.get("intron_start", math.nan)  # relative to gene_start
    intron_end = row.get("intron_end", math.nan)        # relative to gene_start
    
    if not gs or pd.isna(gene_start) or pd.isna(intron_start) or pd.isna(intron_end):
        return "", gene_start, None
    try:
        gene_start = int(gene_start)
        intron_start = int(intron_start)
        intron_end = int(intron_end)
    except (ValueError, TypeError):
        return "", gene_start, None

    # Determine sorting order based on strand.
    strand = row.get("strand", 1)
    if strand in ("-", "-1", -1):
        exons_sorted = sorted(exons, key=lambda e: e["start"], reverse=True)
    else:
        exons_sorted = sorted(exons, key=lambda e: e["start"])
    
    pieces = []
    ir_intron_end_index = None
    for exon in exons_sorted:
        try:
            # Convert exon absolute coordinates to positions relative to gene_start.
            e_start_rel = int(exon["start"] - gene_start)
            e_end_rel = int(exon["end"] - gene_start)
        except (ValueError, TypeError):
            continue
        
        # Extract exon sequence from the gene sequence.
        if 0 <= e_start_rel < e_end_rel <= len(gs):
            exon_seq = gs[e_start_rel:e_end_rel]
        else:
            exon_seq = ""
        pieces.append(exon_seq)
        
        # When the exon ends at the expected intron start, insert the intron sequence.
        if e_end_rel == intron_start:
            pieces.append(intron_seq)
            # Record the IR transcript index at which the intron ends.
            ir_intron_end_index = sum(len(piece) for piece in pieces)
    ir_seq = "".join(pieces)
    return ir_seq, gene_start, ir_intron_end_index

def _translate_dna_modified(dna):
    """
    Translate the DNA sequence beginning at the first ATG.
    Returns a tuple:
      (protein_seq, atg_index)
    where atg_index is the nucleotide offset in 'dna' where translation starts.
    """
    if not dna:
        return "", None
    dna_upper = dna.upper()
    atg_index = dna_upper.find("ATG")
    if atg_index == -1:
        # No ATG found; translate from beginning.
        return str(Seq(dna).translate(to_stop=False)), 0
    trimmed = dna[atg_index:]
    protein_seq = str(Seq(trimmed).translate(to_stop=False))
    return protein_seq, atg_index

def find_first_stop_protein_index(protein_seq):
    """
    Find the index (in amino acids) of the first stop codon in protein_seq.
    Returns None if no stop codon is found.
    """
    stop_pos = protein_seq.find("*")
    if stop_pos == -1:
        return None
    return stop_pos

def find_stop_codon_distance(protein_seq, atg_index, ir_intron_end_index, gene_start, exons, strand):
    """
    Compute two distances (in amino acids) from the first stop codon:
      1. Using the IR transcript:
         The distance from the stop codon to the end of the retained intron.
         (Computed as: (ir_intron_end - (atg_index + stop_codon_nt_offset)) // 3)
      2. Using canonical exon coordinates:
         The distance from the stop codon (projected onto the gene) to the next exon junction.
         (Computed by converting the next exon's absolute start (relative to gene_start) into protein coordinates.)
    
    Parameters:
      protein_seq: The translated protein sequence.
      atg_index: Nucleotide offset in the IR transcript where translation starts.
      ir_intron_end_index: Nucleotide index in the IR transcript where the intron ends.
      gene_start: Absolute gene start (from the Desc field).
      exons: List of exon dictionaries (with absolute 'start' values).
      strand: Strand indicator.
    
    Returns:
      (distance_in_ir, distance_canonical) in amino acids relative to the first stop codon.
    """
    first_stop = find_first_stop_protein_index(protein_seq)
    if first_stop is None:
        return None, None
    # Compute the stop codon nucleotide position in the IR transcript.
    stop_nt_pos = atg_index + first_stop * 3  # offset from beginning of IR transcript
    
    # Distance in the IR transcript:
    if ir_intron_end_index is not None:
        dist_ir = (ir_intron_end_index - stop_nt_pos) // 3
    else:
        dist_ir = None

    # Now compute canonical distance using exons (in absolute coordinates).
    # We project the stop codon back to an approximate absolute position.
    if strand in ("-", "-1", -1):
        stop_abs = gene_start - (first_stop * 3)
        sorted_exons = sorted(exons, key=lambda e: e["start"], reverse=True)
        next_exon = None
        for exon in sorted_exons:
            if exon["start"] < stop_abs:
                next_exon = exon
                break
        if next_exon:
            next_exon_rel = next_exon["start"] - gene_start
            next_exon_prot = next_exon_rel // 3
            canon_dist = next_exon_prot - first_stop
        else:
            canon_dist = None
    else:
        stop_abs = gene_start + (first_stop * 3)
        sorted_exons = sorted(exons, key=lambda e: e["start"])
        next_exon = None
        for exon in sorted_exons:
            if exon["start"] > stop_abs:
                next_exon = exon
                break
        if next_exon:
            next_exon_rel = next_exon["start"] - gene_start
            next_exon_prot = next_exon_rel // 3
            canon_dist = next_exon_prot - first_stop
        else:
            canon_dist = None

    return dist_ir, canon_dist

if __name__ == "__main__":
    # Adjust these file paths and sheet names as needed.
    input_excel = "data/final_output.xlsx"
    output_excel = "data/final_output_updated.xlsx"
    exon_sheet = "Exon Data"
    master_sheet = "Master"
    timepoint_sheets = ["CTRL", "1HR", "2HR", "4HR"]
    main()

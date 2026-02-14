"""
Build a BERDL SQLite database from pipeline output files.

Usage:
    from KBDatalakeApps.build_berdl_db import build_berdl_database
    build_berdl_database("/path/to/pipeline/output", "/path/to/output.db")
"""
import csv
import os
import re
import json
import sqlite3
import hashlib
import glob as globmod

import pandas as pd


def _find_clade_dir(output_dir):
    """Find the pangenome clade directory from user_to_clade.json."""
    clade_file = os.path.join(output_dir, "pangenome", "user_to_clade.json")
    if not os.path.exists(clade_file):
        return None, None
    with open(clade_file) as f:
        mapping = json.load(f)
    # mapping is {genome_name: clade_id}
    genome_name = list(mapping.keys())[0]
    clade_id = mapping[genome_name]
    clade_dir = os.path.join(output_dir, "pangenome", clade_id)
    if os.path.exists(clade_dir):
        return clade_id, clade_dir
    return clade_id, None


def _count_fasta_contigs(fna_path):
    """Count the number of sequences (contigs) in a FASTA file."""
    count = 0
    with open(fna_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def _count_tsv_features(tsv_path):
    """Count features in a genome TSV file (subtract 1 for header)."""
    count = 0
    with open(tsv_path) as f:
        for _ in f:
            count += 1
    return max(count - 1, 0)


def _parse_fasta_sequences(faa_path):
    """Parse a FASTA file into {header_id: sequence} dict."""
    sequences = {}
    current_id = None
    current_seq = []
    with open(faa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)
    return sequences


def _user_genome_id(output_dir):
    """Derive the user genome ID from the kbasedump filename."""
    genome_dir = os.path.join(output_dir, "genome")
    for f in os.listdir(genome_dir):
        if f.endswith("_kbasedump.tsv"):
            return f.replace("_kbasedump.tsv", "")
    return None


def _parse_ec_from_function(func_str):
    """Extract EC numbers from a RAST function string like 'Name (EC x.x.x.x)'."""
    if not func_str or pd.isna(func_str):
        return None
    ecs = re.findall(r'EC[:\s]*([\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+)', str(func_str))
    return ";".join(ecs) if ecs else None


def _parse_gene_name_from_aliases(aliases_str):
    """Extract gene name from aliases like 'protein_id:WP_000831329.1;gene:rpmH;old_locus_tag:...'."""
    if not aliases_str or pd.isna(aliases_str):
        return None
    # Extract just the gene symbol (e.g., "rpmH" from "gene:rpmH")
    match = re.search(r'gene:([^;]+)', str(aliases_str))
    return match.group(1).strip() if match else None


def _parse_ko_annotation(ko_str):
    """Parse KO annotation column: 'KO:K02914:description' → 'K02914'."""
    if not ko_str or pd.isna(ko_str):
        return None
    terms = []
    for part in str(ko_str).split("|"):
        match = re.match(r'KO:(K\d+)', part.strip())
        if match:
            terms.append(match.group(1))
    return ";".join(terms) if terms else None


def _parse_pfam_annotation(pfam_str):
    """Parse PFAM annotation column: 'PFAM:PF:PF00468:desc' → 'PF00468'."""
    if not pfam_str or pd.isna(pfam_str):
        return None
    terms = []
    for part in str(pfam_str).split("|"):
        matches = re.findall(r'(PF\d+)', part)
        terms.extend(matches)
    # Deduplicate while preserving order
    seen = set()
    unique = []
    for t in terms:
        if t not in seen:
            seen.add(t)
            unique.append(t)
    return ";".join(unique) if unique else None


def _parse_sso_annotation(sso_str):
    """Parse SSO annotation column: 'SSO:000004300:desc' → 'SSO:000004300'."""
    if not sso_str or pd.isna(sso_str):
        return None
    terms = []
    for part in str(sso_str).split("|"):
        match = re.match(r'(SSO:\d+)', part.strip())
        if match:
            terms.append(match.group(1))
    return ";".join(terms) if terms else None


def _sha256_hash(seq):
    """Compute SHA-256 hash of a protein sequence."""
    if not seq:
        return ""
    return hashlib.sha256(seq.encode()).hexdigest()


def _parse_bakta_json(bakta_path):
    """Parse a bakta annotation file (Python repr dict format) into a dict keyed by feature ID.

    Returns: {feature_id: {"bakta_function": ..., "cog": ..., "ec": ..., "go": ...,
              "so": ..., "uniref_100": ..., "uniref_90": ..., "uniref_50": ..., "gene": ...}}
    """
    import ast
    if not os.path.exists(bakta_path):
        return {}
    with open(bakta_path) as f:
        raw = f.read()
    # Strip Docker log timestamp corruption
    raw = re.sub(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z?\s*', '', raw)
    try:
        data = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        try:
            data = json.loads(raw)
        except json.JSONDecodeError:
            print(f"WARNING: Could not parse bakta file: {bakta_path}")
            return {}

    result = {}
    for feat in data.get("features", []):
        fid = feat.get("id") or feat.get("locus", "")
        if not fid:
            continue

        psc = feat.get("psc", {}) or {}
        ips = feat.get("ips", {}) or {}
        ups = feat.get("ups", {}) or {}
        pscc = feat.get("pscc", {}) or {}

        # COG
        cog_id = psc.get("cog_id", "")

        # EC from psc.ec_ids
        ec_ids = psc.get("ec_ids", []) or []
        ec_str = ";".join(ec_ids) if ec_ids else ""

        # GO from psc.go_ids
        go_ids = psc.get("go_ids", []) or []
        go_str = ";".join(sorted(go_ids)) if go_ids else ""

        # SO from db_xrefs
        db_xrefs = feat.get("db_xrefs", []) or []
        so_terms = sorted(x for x in db_xrefs if x.startswith("SO:"))
        so_str = ";".join(so_terms) if so_terms else ""

        # UniRef from various sub-objects
        uniref100 = ups.get("uniref100_id", "") or ips.get("uniref100_id", "") or ""
        uniref90 = ips.get("uniref90_id", "") or psc.get("uniref90_id", "") or ""
        uniref50 = psc.get("uniref50_id", "") or pscc.get("uniref50_id", "") or ""

        result[fid] = {
            "bakta_function": feat.get("product", ""),
            "cog": cog_id,
            "ec": ec_str,
            "go": go_str,
            "so": so_str,
            "uniref_100": uniref100,
            "uniref_90": uniref90,
            "uniref_50": uniref50,
            "gene": feat.get("gene", ""),
        }
    return result


def _parse_kofamscan_tsv(kofam_path):
    """Parse KOfamscan TSV (id, functions) into {feature_id: ko_id}."""
    if not os.path.exists(kofam_path):
        return {}
    result = {}
    with open(kofam_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fid = row.get("id", "").strip()
            ko = row.get("functions", "").strip()
            if fid and ko:
                result[fid] = ko
    return result


def _parse_psortb_tsv(psort_path):
    """Parse PSORTb TSV into {feature_id: localization_string}."""
    if not os.path.exists(psort_path):
        return {}
    result = {}
    with open(psort_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fid = row.get("SeqID", "").strip()
            loc = row.get("CMSVM-_Localization", "").strip()
            if fid and loc and loc.lower() != "unknown":
                result[fid] = loc
    return result


# ---------------------------------------------------------------------------
# Table 1: genome
# ---------------------------------------------------------------------------
def build_genome_table(output_dir, conn):
    """
    Build the genome table with metadata for user + pangenome member genomes.

    Columns: id, gtdb_taxonomy, ncbi_taxonomy, n_contigs, n_features
    """
    rows = []
    user_gid = _user_genome_id(output_dir)
    clade_id, clade_dir = _find_clade_dir(output_dir)

    # --- User genome ---
    user_row = {"id": user_gid, "gtdb_taxonomy": None, "ncbi_taxonomy": None,
                "n_contigs": 0, "n_features": 0}
    # Count contigs from assembly
    assembly_dir = os.path.join(output_dir, "assembly")
    if os.path.exists(assembly_dir):
        for fna in os.listdir(assembly_dir):
            if fna.endswith(".fna"):
                user_row["n_contigs"] = _count_fasta_contigs(
                    os.path.join(assembly_dir, fna))
                break
    # Count features from kbasedump
    genome_dir = os.path.join(output_dir, "genome")
    for f in os.listdir(genome_dir):
        if f.endswith("_kbasedump.tsv"):
            user_row["n_features"] = _count_tsv_features(
                os.path.join(genome_dir, f))
            break
    rows.append(user_row)

    # --- Pangenome member genomes ---
    if clade_dir:
        members_file = os.path.join(clade_dir, "members.tsv")
        if os.path.exists(members_file):
            members_df = pd.read_csv(members_file, sep="\t")
            pan_assembly_dir = os.path.join(clade_dir, "assembly")
            pan_genome_dir = os.path.join(clade_dir, "genome")
            for _, member in members_df.iterrows():
                gid = member["genome_id"]
                row = {
                    "id": gid,
                    "gtdb_taxonomy": member.get("gtdb_taxonomy_id", None),
                    "ncbi_taxonomy": None,  # MISSING: not in input files
                    "n_contigs": 0,
                    "n_features": 0,
                }
                # Count contigs
                fna = os.path.join(pan_assembly_dir, f"{gid}.fna")
                if os.path.exists(fna):
                    row["n_contigs"] = _count_fasta_contigs(fna)
                # Count features
                tsv = os.path.join(pan_genome_dir, f"{gid}.tsv")
                if os.path.exists(tsv):
                    row["n_features"] = _count_tsv_features(tsv)
                rows.append(row)

    df = pd.DataFrame(rows)
    # Ensure column order matches target schema
    for col in ["id", "gtdb_taxonomy", "ncbi_taxonomy", "n_contigs", "n_features"]:
        if col not in df.columns:
            df[col] = None

    df = df[["id", "gtdb_taxonomy", "ncbi_taxonomy", "n_contigs", "n_features"]]
    df.to_sql("genome", conn, if_exists="replace", index=False,
              dtype={
                  "id": "VARCHAR(255) NOT NULL",
                  "gtdb_taxonomy": "VARCHAR(1000)",
                  "ncbi_taxonomy": "VARCHAR(1000)",
                  "n_contigs": "INTEGER",
                  "n_features": "INTEGER",
              })
    # Add primary key via ALTER is not supported in SQLite; use CREATE TABLE approach
    # The schema is created by to_sql; we rely on the VARCHAR(255) NOT NULL for id
    print(f"Built 'genome' table: {len(df)} rows")
    return df


# ---------------------------------------------------------------------------
# Table 2: genome_ani
# ---------------------------------------------------------------------------
def build_genome_ani_table(output_dir, conn):
    """
    Build the genome_ani table from ANI output files.

    Columns: genome1, genome2, ani, af1, af2, kind
    """
    ani_dir = os.path.join(output_dir, "ani")
    if not os.path.exists(ani_dir):
        print("No ani/ directory found, skipping genome_ani table")
        return pd.DataFrame()

    user_gid = _user_genome_id(output_dir)

    def _extract_genome_id_from_path(filepath_str):
        """Extract genome ID from a file path like '.../GCF_000368685.1_.../file.fna'."""
        if not filepath_str or pd.isna(filepath_str):
            return ""
        s = str(filepath_str)
        # Check for user genome in path
        if "user_" in s:
            match = re.search(r'(user_[^/]+?)\.fna', s)
            if match:
                return match.group(1)
        # Try to extract GCF/GCA ID from path
        match = re.search(r'(GC[AF]_\d+\.\d+)', s)
        if match:
            return "RS_" + match.group(1) if "GCF" in match.group(1) else "GB_" + match.group(1)
        return s

    rows = []
    for ani_file in os.listdir(ani_dir):
        if not ani_file.endswith(".out"):
            continue
        kind = ani_file.replace("_fast.out", "").replace(".out", "")
        filepath = os.path.join(ani_dir, ani_file)
        try:
            df = pd.read_csv(filepath, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if len(df) == 0:
            continue
        for _, row in df.iterrows():
            # Use file paths to get proper genome IDs, falling back to name columns
            genome1 = _extract_genome_id_from_path(row.get("Query_file", ""))
            genome2 = _extract_genome_id_from_path(row.get("Ref_file", ""))
            if not genome1:
                genome1 = row.get("Query_name", row.get("query_name", ""))
            if not genome2:
                genome2 = row.get("Ref_name", row.get("ref_name", ""))
            rows.append({
                "genome1": genome1,
                "genome2": genome2,
                "ani": row.get("ANI", row.get("ani", 0.0)),
                "af1": row.get("Align_fraction_ref", row.get("align_fraction_ref", 0.0)),
                "af2": row.get("Align_fraction_query", row.get("align_fraction_query", 0.0)),
                "kind": kind,
            })

    if not rows:
        print("No ANI data found, skipping genome_ani table")
        return pd.DataFrame()

    result_df = pd.DataFrame(rows)
    result_df.to_sql("genome_ani", conn, if_exists="replace", index=False,
                     dtype={
                         "genome1": "VARCHAR(255) NOT NULL",
                         "genome2": "VARCHAR(255) NOT NULL",
                         "ani": "FLOAT NOT NULL",
                         "af1": "FLOAT NOT NULL",
                         "af2": "FLOAT NOT NULL",
                         "kind": "VARCHAR(255) NOT NULL",
                     })
    print(f"Built 'genome_ani' table: {len(result_df)} rows")
    return result_df


# ---------------------------------------------------------------------------
# Table 3: genome_features
# ---------------------------------------------------------------------------
def build_genome_features_table(output_dir, conn):
    """
    Build the genome_features table for the user genome.

    Main source: genome/*_kbasedump.tsv
    Additional sources: pangenome cluster parquet, gene_reaction_data.tsv
    """
    user_gid = _user_genome_id(output_dir)
    genome_dir = os.path.join(output_dir, "genome")

    # Load kbasedump TSV
    kbasedump_path = os.path.join(genome_dir, f"{user_gid}_kbasedump.tsv")
    if not os.path.exists(kbasedump_path):
        print(f"kbasedump file not found: {kbasedump_path}")
        return pd.DataFrame()

    kdf = pd.read_csv(kbasedump_path, sep="\t")

    # Build features dataframe
    features = pd.DataFrame()
    features["genome_id"] = [user_gid] * len(kdf)
    features["contig_id"] = kdf["contig"].values
    features["feature_id"] = kdf["gene_id"].values
    features["start"] = kdf["start"].values
    features["end"] = kdf["end"].values
    features["length"] = abs(kdf["end"] - kdf["start"]).values
    features["strand"] = kdf["strand"].values
    features["sequence"] = kdf["protein_translation"].values if "protein_translation" in kdf.columns else None
    features["sequence_hash"] = features["sequence"].apply(
        lambda x: hashlib.md5(str(x).encode()).hexdigest() if pd.notna(x) and x else "")

    # RAST function
    features["rast_function"] = kdf["functions"].values

    # Parse annotations from kbasedump columns
    ko_series = kdf["Annotation:KO"] if "Annotation:KO" in kdf.columns else pd.Series([None]*len(kdf))
    pfam_series = kdf["Annotation:PFAM"] if "Annotation:PFAM" in kdf.columns else pd.Series([None]*len(kdf))
    sso_series = kdf["Annotation:SSO"] if "Annotation:SSO" in kdf.columns else pd.Series([None]*len(kdf))
    features["ko"] = ko_series.apply(_parse_ko_annotation).values
    features["pfam"] = pfam_series.apply(_parse_pfam_annotation).values
    features["so"] = sso_series.apply(_parse_sso_annotation).values

    # Extract EC from RAST function strings
    features["ec"] = kdf["functions"].apply(_parse_ec_from_function).values

    # Extract gene names from aliases
    alias_series = kdf["aliases"] if "aliases" in kdf.columns else pd.Series([None]*len(kdf))
    features["gene_names"] = alias_series.apply(_parse_gene_name_from_aliases).values

    # --- Load Bakta annotations if available ---
    bakta_data = {}
    for ext in [".json", ".tsv"]:
        bakta_path = os.path.join(genome_dir, f"{user_gid}_bakta{ext}")
        if os.path.exists(bakta_path):
            bakta_data = _parse_bakta_json(bakta_path)
            if bakta_data:
                print(f"Loaded {len(bakta_data)} Bakta annotations from {os.path.basename(bakta_path)}")
            break

    if bakta_data:
        features["bakta_function"] = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("bakta_function") or None)
        features["cog"] = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("cog") or None)
        features["go"] = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("go") or None)
        features["uniref_100"] = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("uniref_100") or None)
        features["uniref_90"] = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("uniref_90") or None)
        features["uniref_50"] = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("uniref_50") or None)
        # Merge bakta EC with RAST-extracted EC
        bakta_ec = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("ec") or None)
        features["ec"] = features.apply(
            lambda r: ";".join(filter(None, set((str(r["ec"]) if pd.notna(r["ec"]) else "").split(";") +
                                                (str(bakta_ec[r.name]) if pd.notna(bakta_ec[r.name]) else "").split(";")))) or None,
            axis=1)
        # Merge bakta SO with kbasedump SO
        bakta_so = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("so") or None)
        features["so"] = features.apply(
            lambda r: ";".join(sorted(set(filter(None,
                (str(r["so"]) if pd.notna(r["so"]) else "").split(";") +
                (str(bakta_so[r.name]) if pd.notna(bakta_so[r.name]) else "").split(";"))))) or None,
            axis=1)
        # Merge bakta gene names with parsed aliases
        bakta_gene = features["feature_id"].map(
            lambda fid: bakta_data.get(fid, {}).get("gene") or None)
        features["gene_names"] = features.apply(
            lambda r: r["gene_names"] if pd.notna(r["gene_names"]) else bakta_gene[r.name],
            axis=1)
    else:
        features["bakta_function"] = None
        features["cog"] = None
        features["go"] = None
        features["uniref_100"] = None
        features["uniref_90"] = None
        features["uniref_50"] = None

    # --- Load KOfamscan annotations if available ---
    kofam_path = os.path.join(genome_dir, f"{user_gid}_KOfamscan.tsv")
    kofam_data = _parse_kofamscan_tsv(kofam_path)
    if kofam_data:
        print(f"Loaded {len(kofam_data)} KOfamscan annotations")
        # Merge KOfamscan KO with kbasedump KO
        kofam_ko = features["feature_id"].map(lambda fid: kofam_data.get(fid))
        features["ko"] = features.apply(
            lambda r: ";".join(sorted(set(filter(None,
                (str(r["ko"]) if pd.notna(r["ko"]) else "").split(";") +
                (str(kofam_ko[r.name]) if pd.notna(kofam_ko[r.name]) else "").split(";"))))) or None,
            axis=1)

    # --- Load PSORTb annotations if available ---
    psort_path = os.path.join(genome_dir, f"{user_gid}_PSORT.tsv")
    psortb_data = _parse_psortb_tsv(psort_path)
    if psortb_data:
        print(f"Loaded {len(psortb_data)} PSORTb annotations")
        features["psortb"] = features["feature_id"].map(lambda fid: psortb_data.get(fid))
    else:
        features["psortb"] = None

    features["rast_consistency"] = None  # Not available
    features["other_rast_annotations"] = None  # Not available

    # --- Pangenome cluster assignment ---
    # The user genome is NOT in the pangenome parquet directly.
    # We need to: hash user proteins → look up in mmseqs2 cluster.tsv → match to cluster IDs
    clade_id, clade_dir = _find_clade_dir(output_dir)
    features["pangenome_cluster_id"] = None
    features["pangenome_is_core"] = None

    if clade_dir:
        pangenome_parquet = os.path.join(clade_dir, "pangenome_cluster_with_mmseqs.parquet")
        cluster_tsv = os.path.join(clade_dir, "master_mmseqs2", "master_faa_cluster.tsv")

        if os.path.exists(pangenome_parquet) and os.path.exists(cluster_tsv):
            pan_df = pd.read_parquet(pangenome_parquet)

            # Build rep_hash → cluster_id mapping from pangenome parquet
            rep_to_cluster = pan_df.drop_duplicates("mmseqs_rep_hash").set_index(
                "mmseqs_rep_hash")["cluster_id"].to_dict()

            # Compute core clusters (present in all pangenome genomes)
            n_pan_genomes = pan_df["genome_id"].nunique()
            cluster_genome_counts = pan_df.groupby("cluster_id")["genome_id"].nunique()
            core_clusters = set(cluster_genome_counts[cluster_genome_counts == n_pan_genomes].index)

            # Load mmseqs2 cluster assignments: member_hash → rep_hash
            cluster_map = {}
            with open(cluster_tsv) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        rep_hash, member_hash = parts
                        cluster_map[member_hash] = rep_hash

            # Hash user genome protein sequences and look up clusters
            faa_path = os.path.join(genome_dir, f"{user_gid}.faa")
            if os.path.exists(faa_path):
                user_seqs = _parse_fasta_sequences(faa_path)
                # Map feature_id → protein hash → rep hash → cluster_id
                feature_to_cluster = {}
                for feat_id, seq in user_seqs.items():
                    protein_hash = _sha256_hash(seq)
                    rep_hash = cluster_map.get(protein_hash)
                    if rep_hash:
                        cluster_id = rep_to_cluster.get(rep_hash)
                        if cluster_id:
                            feature_to_cluster[feat_id] = cluster_id

                features["pangenome_cluster_id"] = features["feature_id"].map(feature_to_cluster)
                features["pangenome_is_core"] = features["pangenome_cluster_id"].apply(
                    lambda x: x in core_clusters if pd.notna(x) else None)

                assigned = features["pangenome_cluster_id"].notna().sum()
                total = len(features)
                print(f"Pangenome cluster assignment: {assigned}/{total} features assigned")

    # --- Reactions and model data from gene_reaction_data.tsv ---
    gene_rxn_path = os.path.join(output_dir, "models", "gene_reaction_data.tsv")
    if os.path.exists(gene_rxn_path):
        gene_rxn_df = pd.read_csv(gene_rxn_path, sep="\t")
        user_rxn = gene_rxn_df[gene_rxn_df["genome_id"] == user_gid]
        if len(user_rxn) > 0:
            rxn_idx = user_rxn.set_index("gene_id")
            features["reactions"] = features["feature_id"].map(rxn_idx["reaction"].to_dict())
            features["rich_media_flux"] = features["feature_id"].map(rxn_idx["rich_media_flux"].to_dict())
            features["rich_media_class"] = features["feature_id"].map(rxn_idx["rich_media_class"].to_dict())
            features["minimal_media_flux"] = features["feature_id"].map(rxn_idx["minimal_media_flux"].to_dict())
            features["minimal_media_class"] = features["feature_id"].map(rxn_idx["minimal_media_class"].to_dict())
        else:
            features["reactions"] = None
            features["rich_media_flux"] = None
            features["rich_media_class"] = None
            features["minimal_media_flux"] = None
            features["minimal_media_class"] = None
    else:
        features["reactions"] = None
        features["rich_media_flux"] = None
        features["rich_media_class"] = None
        features["minimal_media_flux"] = None
        features["minimal_media_class"] = None

    # Reorder columns to match target schema
    col_order = [
        "genome_id", "contig_id", "feature_id", "length", "start", "end",
        "strand", "sequence", "sequence_hash", "bakta_function", "rast_function",
        "cog", "ec", "gene_names", "go", "ko", "pfam", "so",
        "uniref_100", "uniref_90", "uniref_50",
        "pangenome_cluster_id", "pangenome_is_core",
        "psortb", "reactions", "rich_media_flux", "rich_media_class",
        "minimal_media_flux", "minimal_media_class",
        "rast_consistency", "other_rast_annotations",
    ]
    features = features[col_order]

    # Write to SQLite with auto-increment id
    features.to_sql("genome_features", conn, if_exists="replace", index=True,
                    index_label="id",
                    dtype={
                        "id": "INTEGER NOT NULL",
                        "genome_id": "VARCHAR(255) NOT NULL",
                        "contig_id": "VARCHAR(255) NOT NULL",
                        "feature_id": "VARCHAR(255) NOT NULL",
                        "length": "INTEGER NOT NULL",
                        "start": "INTEGER",
                        "end": "INTEGER",
                        "strand": "VARCHAR(1)",
                        "sequence": "TEXT",
                        "sequence_hash": "VARCHAR(255)",
                        "bakta_function": "VARCHAR(255)",
                        "rast_function": "VARCHAR(255)",
                        "cog": "VARCHAR(255)",
                        "ec": "VARCHAR(255)",
                        "gene_names": "VARCHAR(255)",
                        "go": "VARCHAR(1000)",
                        "ko": "VARCHAR(1000)",
                        "pfam": "VARCHAR(1000)",
                        "so": "VARCHAR(1000)",
                        "uniref_100": "VARCHAR(255)",
                        "uniref_90": "VARCHAR(255)",
                        "uniref_50": "VARCHAR(255)",
                        "pangenome_cluster_id": "VARCHAR(255)",
                        "pangenome_is_core": "BOOLEAN",
                        "psortb": "VARCHAR(255)",
                        "reactions": "TEXT",
                        "rich_media_flux": "REAL",
                        "rich_media_class": "VARCHAR(255)",
                        "minimal_media_flux": "REAL",
                        "minimal_media_class": "VARCHAR(255)",
                        "rast_consistency": "REAL",
                        "other_rast_annotations": "TEXT",
                    })

    # Add unique constraint
    try:
        conn.execute(
            "CREATE UNIQUE INDEX IF NOT EXISTS idx_genome_features_unique "
            "ON genome_features(genome_id, contig_id, feature_id)")
    except sqlite3.OperationalError:
        pass

    print(f"Built 'genome_features' table: {len(features)} rows")
    return features


# ---------------------------------------------------------------------------
# Table 4: missing_functions
# ---------------------------------------------------------------------------
def build_missing_functions_table(output_dir, conn):
    """
    Build the missing_functions table from gapfilling data.

    Sources: genome_reactions.tsv, model *_data.json, genome_phenotypes.tsv,
             and pangenome comparison.
    """
    user_gid = _user_genome_id(output_dir)
    all_reactions = {}  # reaction_id → {column: value}

    # --- Gapfilled reactions from genome_reactions.tsv ---
    rxn_path = os.path.join(output_dir, "models", "genome_reactions.tsv")
    if os.path.exists(rxn_path):
        rxn_df = pd.read_csv(rxn_path, sep="\t")
        user_rxn = rxn_df[rxn_df["genome_id"] == user_gid]
        gapfilled = user_rxn[user_rxn["gapfilling_status"] != "none"]
        for _, row in gapfilled.iterrows():
            rxn_id = row["reaction_id"]
            status = row["gapfilling_status"]
            if rxn_id not in all_reactions:
                all_reactions[rxn_id] = {
                    "Reaction": rxn_id, "RAST_function": None,
                    "RichGapfill": 0, "MinimalGapfill": 0,
                    "PhenotypeGapfill": 0, "ModuleGapfill": 0, "Pangenome": 0,
                }
            if "rich" in status.lower():
                all_reactions[rxn_id]["RichGapfill"] = 1
            if "minimal" in status.lower():
                all_reactions[rxn_id]["MinimalGapfill"] = 1
            if "core" in status.lower():
                all_reactions[rxn_id]["MinimalGapfill"] = 1

    # --- Gapfilled reactions from model *_data.json files ---
    models_dir = os.path.join(output_dir, "models")
    user_data_files = globmod.glob(os.path.join(models_dir, f"{user_gid}_data.json"))
    for data_file in user_data_files:
        with open(data_file) as f:
            model_data = json.load(f)
        gf = model_data.get("gapfilled_reactions", {})
        for category, rxn_list in gf.items():
            for rxn_id in rxn_list:
                if rxn_id not in all_reactions:
                    all_reactions[rxn_id] = {
                        "Reaction": rxn_id, "RAST_function": None,
                        "RichGapfill": 0, "MinimalGapfill": 0,
                        "PhenotypeGapfill": 0, "ModuleGapfill": 0, "Pangenome": 0,
                    }
                if "rich" in category.lower():
                    all_reactions[rxn_id]["RichGapfill"] = 1
                if "minimal" in category.lower():
                    all_reactions[rxn_id]["MinimalGapfill"] = 1

    # --- Phenotype gapfilled reactions from genome_phenotypes.tsv ---
    pheno_path = os.path.join(output_dir, "phenotype_tables", "genome_phenotypes.tsv")
    if os.path.exists(pheno_path):
        pheno_df = pd.read_csv(pheno_path, sep="\t")
        user_pheno = pheno_df[pheno_df["genome_id"] == user_gid]
        for _, row in user_pheno.iterrows():
            gf_rxns = row.get("gapfilled_reactions", "")
            if pd.notna(gf_rxns) and gf_rxns:
                for rxn_id in str(gf_rxns).split(";"):
                    rxn_id = rxn_id.strip()
                    if not rxn_id:
                        continue
                    if rxn_id not in all_reactions:
                        all_reactions[rxn_id] = {
                            "Reaction": rxn_id, "RAST_function": None,
                            "RichGapfill": 0, "MinimalGapfill": 0,
                            "PhenotypeGapfill": 0, "ModuleGapfill": 0, "Pangenome": 0,
                        }
                    all_reactions[rxn_id]["PhenotypeGapfill"] = 1

    # --- Pangenome-derived missing reactions ---
    # Reactions present in pangenome member models but not in user model
    if os.path.exists(rxn_path):
        rxn_df = pd.read_csv(rxn_path, sep="\t")
        user_rxns = set(rxn_df[rxn_df["genome_id"] == user_gid]["reaction_id"])
        pan_rxns = set(rxn_df[rxn_df["genome_id"] != user_gid]["reaction_id"])
        missing_in_user = pan_rxns - user_rxns
        for rxn_id in missing_in_user:
            if rxn_id not in all_reactions:
                all_reactions[rxn_id] = {
                    "Reaction": rxn_id, "RAST_function": None,
                    "RichGapfill": 0, "MinimalGapfill": 0,
                    "PhenotypeGapfill": 0, "ModuleGapfill": 0, "Pangenome": 1,
                }
            else:
                all_reactions[rxn_id]["Pangenome"] = 1

    # --- Try to map reactions to RAST functions ---
    # Use genome_reactions.tsv which has both reaction_id and genes
    if os.path.exists(rxn_path) and all_reactions:
        rxn_df = pd.read_csv(rxn_path, sep="\t")
        # Build reaction → equation_names mapping
        rxn_to_name = rxn_df.drop_duplicates("reaction_id").set_index(
            "reaction_id")["equation_names"].to_dict()
        for rxn_id, entry in all_reactions.items():
            if entry["RAST_function"] is None:
                entry["RAST_function"] = rxn_to_name.get(rxn_id, None)

    if not all_reactions:
        print("No missing functions found, creating empty missing_functions table")
        empty_df = pd.DataFrame(columns=[
            "Reaction", "RAST_function", "RichGapfill", "MinimalGapfill",
            "PhenotypeGapfill", "ModuleGapfill", "Pangenome"])
        empty_df.to_sql("missing_functions", conn, if_exists="replace", index=False)
        return empty_df

    result_df = pd.DataFrame(list(all_reactions.values()))
    result_df = result_df[["Reaction", "RAST_function", "RichGapfill",
                           "MinimalGapfill", "PhenotypeGapfill",
                           "ModuleGapfill", "Pangenome"]]
    result_df.to_sql("missing_functions", conn, if_exists="replace", index=False,
                     dtype={
                         "Reaction": "TEXT",
                         "RAST_function": "TEXT",
                         "RichGapfill": "INTEGER",
                         "MinimalGapfill": "INTEGER",
                         "PhenotypeGapfill": "INTEGER",
                         "ModuleGapfill": "INTEGER",
                         "Pangenome": "INTEGER",
                     })
    print(f"Built 'missing_functions' table: {len(result_df)} rows")
    return result_df


# ---------------------------------------------------------------------------
# Table 5: pan_genome_features
# ---------------------------------------------------------------------------
def build_pan_genome_features_table(output_dir, conn):
    """
    Build the pan_genome_features table from pangenome member genome files.

    Sources: pangenome/{clade}/genome/*.tsv, *.faa, pangenome_cluster_with_mmseqs.parquet
    """
    clade_id, clade_dir = _find_clade_dir(output_dir)
    if not clade_dir:
        print("No pangenome clade directory found, skipping pan_genome_features table")
        return pd.DataFrame()

    pan_genome_dir = os.path.join(clade_dir, "genome")
    pangenome_parquet = os.path.join(clade_dir, "pangenome_cluster_with_mmseqs.parquet")

    if not os.path.exists(pan_genome_dir):
        print("No pangenome genome directory found")
        return pd.DataFrame()

    # Load cluster data
    cluster_df = None
    core_clusters = set()
    if os.path.exists(pangenome_parquet):
        cluster_df = pd.read_parquet(pangenome_parquet)
        n_genomes = cluster_df["genome_id"].nunique()
        cluster_genome_counts = cluster_df.groupby("cluster_id")["genome_id"].nunique()
        core_clusters = set(cluster_genome_counts[cluster_genome_counts == n_genomes].index)
        # Build (genome_id, feature_id) → cluster_id map
        cluster_map = {}
        for _, row in cluster_df.iterrows():
            cluster_map[(row["genome_id"], row["feature_id"])] = row["cluster_id"]

    all_features = []

    for tsv_file in sorted(os.listdir(pan_genome_dir)):
        if not tsv_file.endswith(".tsv") or tsv_file.endswith("_genome_data.tsv"):
            continue
        genome_id = tsv_file.replace(".tsv", "")
        tsv_path = os.path.join(pan_genome_dir, tsv_file)
        faa_path = os.path.join(pan_genome_dir, f"{genome_id}.faa")
        enriched_path = os.path.join(pan_genome_dir, f"{genome_id}_genome_data.tsv")

        # Check if enriched TSV exists (from enrich_pangenome_features.py)
        enriched = {}
        if os.path.exists(enriched_path):
            edf = pd.read_csv(enriched_path, sep="\t", dtype=str).fillna("")
            for _, erow in edf.iterrows():
                enriched[erow["feature_id"]] = erow.to_dict()
            print(f"  Loaded {len(enriched)} enriched features for {genome_id}")

        # Read feature annotations from base RAST TSV
        gdf = pd.read_csv(tsv_path, sep="\t")

        # Read protein sequences
        sequences = {}
        if os.path.exists(faa_path):
            sequences = _parse_fasta_sequences(faa_path)

        for _, row in gdf.iterrows():
            feature_id = row["id"]
            rast_func = row.get("functions", None)

            # Check enriched data for this feature
            edata = enriched.get(feature_id, {})

            # Use enriched contig_id if available, else parse from feature_id
            contig_id = edata.get("contig_id", "")
            if not contig_id:
                parts = feature_id.rsplit("_", 1)
                contig_id = parts[0] if len(parts) == 2 and parts[1].isdigit() else feature_id

            # Sequence: prefer enriched (has full sequence), fallback to FAA
            seq = edata.get("sequence", "") or sequences.get(feature_id, None)
            seq_hash = hashlib.md5(seq.encode()).hexdigest() if seq else ""
            length = int(edata["length"]) if edata.get("length") else (len(seq) * 3 if seq else 0)

            # Coordinates from enriched data
            start = int(edata["start"]) if edata.get("start") else None
            end = int(edata["end"]) if edata.get("end") else None

            cluster_id = None
            is_core = None
            if cluster_df is not None:
                cluster_id = cluster_map.get((genome_id, feature_id))
                if cluster_id:
                    is_core = cluster_id in core_clusters

            # Use enriched RAST function if the base one is empty
            if pd.isna(rast_func) or not rast_func:
                rast_func = edata.get("rast_function", None) or None

            ec = edata.get("ec", "") or _parse_ec_from_function(rast_func)

            all_features.append({
                "genome_id": genome_id,
                "contig_id": contig_id,
                "feature_id": feature_id,
                "length": length,
                "start": start,
                "end": end,
                "strand": None,  # Still not available
                "sequence": seq,
                "sequence_hash": seq_hash,
                "cluster_id": cluster_id,
                "is_core": is_core,
                "bakta_function": edata.get("bakta_function") or None,
                "rast_function": rast_func,
                "gene_names": None,  # Still not available for pangenome members
                "cog": edata.get("cog") or None,
                "ec": ec or None,
                "ko": edata.get("ko") or None,
                "pfam": edata.get("pfam") or None,
                "go": edata.get("go") or None,
                "so": edata.get("so") or None,
                "uniref_100": edata.get("uniref_100") or None,
                "uniref_90": edata.get("uniref_90") or None,
                "uniref_50": edata.get("uniref_50") or None,
            })

    if not all_features:
        print("No pangenome features found")
        return pd.DataFrame()

    result_df = pd.DataFrame(all_features)

    # --- Merge model data from gene_reaction_data.tsv ---
    gene_rxn_path = os.path.join(output_dir, "models", "gene_reaction_data.tsv")
    if os.path.exists(gene_rxn_path):
        gene_rxn_df = pd.read_csv(gene_rxn_path, sep="\t")
        rxn_idx = gene_rxn_df.set_index(["genome_id", "gene_id"])
        keys = list(zip(result_df["genome_id"], result_df["feature_id"]))
        for col in ["reaction", "rich_media_flux", "rich_media_class",
                     "minimal_media_flux", "minimal_media_class"]:
            result_df[col] = [
                rxn_idx.loc[k, col] if k in rxn_idx.index else None
                for k in keys
            ]
    else:
        for col in ["reaction", "rich_media_flux", "rich_media_class",
                     "minimal_media_flux", "minimal_media_class"]:
            result_df[col] = None

    col_order = [
        "genome_id", "contig_id", "feature_id", "length", "start", "end",
        "strand", "sequence", "sequence_hash", "cluster_id", "is_core",
        "bakta_function", "rast_function", "gene_names",
        "cog", "ec", "ko", "pfam", "go", "so",
        "uniref_100", "uniref_90", "uniref_50",
        "reaction", "rich_media_flux", "rich_media_class",
        "minimal_media_flux", "minimal_media_class",
    ]
    result_df = result_df[col_order]

    result_df.to_sql("pan_genome_features", conn, if_exists="replace", index=True,
                     index_label="id",
                     dtype={
                         "id": "INTEGER NOT NULL",
                         "genome_id": "VARCHAR(255) NOT NULL",
                         "contig_id": "VARCHAR(255) NOT NULL",
                         "feature_id": "VARCHAR(255) NOT NULL",
                         "length": "INTEGER NOT NULL",
                         "start": "INTEGER",
                         "end": "INTEGER",
                         "strand": "VARCHAR(1)",
                         "sequence": "TEXT",
                         "sequence_hash": "VARCHAR(255)",
                         "cluster_id": "VARCHAR(255)",
                         "is_core": "BOOLEAN",
                         "bakta_function": "VARCHAR(255)",
                         "rast_function": "VARCHAR(255)",
                         "gene_names": "VARCHAR(255)",
                         "cog": "VARCHAR(1000)",
                         "ec": "VARCHAR(1000)",
                         "ko": "VARCHAR(1000)",
                         "pfam": "VARCHAR(1000)",
                         "go": "VARCHAR(1000)",
                         "so": "VARCHAR(1000)",
                         "uniref_100": "VARCHAR(255)",
                         "uniref_90": "VARCHAR(255)",
                         "uniref_50": "VARCHAR(255)",
                         "reaction": "TEXT",
                         "rich_media_flux": "REAL",
                         "rich_media_class": "VARCHAR(255)",
                         "minimal_media_flux": "REAL",
                         "minimal_media_class": "VARCHAR(255)",
                     })

    try:
        conn.execute(
            "CREATE UNIQUE INDEX IF NOT EXISTS idx_pan_genome_features_unique "
            "ON pan_genome_features(genome_id, contig_id, feature_id)")
    except sqlite3.OperationalError:
        pass

    print(f"Built 'pan_genome_features' table: {len(result_df)} rows")
    return result_df


# ---------------------------------------------------------------------------
# Table 6: gene_reaction_data
# ---------------------------------------------------------------------------
def build_gene_reaction_data_table(output_dir, conn):
    """
    Load gene_reaction_data.tsv as a standalone table in the database.

    Columns: genome_id, gene_id, reaction, rich_media_flux, rich_media_class,
             minimal_media_flux, minimal_media_class
    """
    gene_rxn_path = os.path.join(output_dir, "models", "gene_reaction_data.tsv")
    if not os.path.exists(gene_rxn_path):
        print("No gene_reaction_data.tsv found, skipping gene_reaction_data table")
        return pd.DataFrame()

    df = pd.read_csv(gene_rxn_path, sep="\t")
    df.to_sql("gene_reaction_data", conn, if_exists="replace", index=False,
              dtype={
                  "genome_id": "VARCHAR(255) NOT NULL",
                  "gene_id": "VARCHAR(255) NOT NULL",
                  "reaction": "TEXT",
                  "rich_media_flux": "REAL",
                  "rich_media_class": "VARCHAR(255)",
                  "minimal_media_flux": "REAL",
                  "minimal_media_class": "VARCHAR(255)",
              })
    print(f"Built 'gene_reaction_data' table: {len(df)} rows")
    return df


# ---------------------------------------------------------------------------
# Table 7: genome_reactions
# ---------------------------------------------------------------------------
def build_genome_reactions_table(output_dir, conn):
    """Load genome_reactions.tsv as the genome_reactions table."""
    tsv_path = os.path.join(output_dir, "models", "genome_reactions.tsv")
    if not os.path.exists(tsv_path):
        print("No genome_reactions.tsv found, skipping genome_reactions table")
        return pd.DataFrame()

    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    df.to_sql("genome_reactions", conn, if_exists="replace", index=False)
    print(f"Built 'genome_reactions' table: {len(df)} rows")
    return df


# ---------------------------------------------------------------------------
# Table 8: growth_phenotype_summary
# ---------------------------------------------------------------------------
def build_growth_phenotype_summary_table(output_dir, conn):
    """Load model_performance.tsv as the growth_phenotype_summary table."""
    tsv_path = os.path.join(output_dir, "phenotypes", "model_performance.tsv")
    if not os.path.exists(tsv_path):
        print("No model_performance.tsv found, skipping growth_phenotype_summary table")
        return pd.DataFrame()

    df = pd.read_csv(tsv_path, sep="\t")
    df.to_sql("growth_phenotype_summary", conn, if_exists="replace", index=False)
    print(f"Built 'growth_phenotype_summary' table: {len(df)} rows")
    return df


# ---------------------------------------------------------------------------
# Table 8: growth_phenotypes_detailed
# ---------------------------------------------------------------------------
def build_growth_phenotypes_detailed_table(output_dir, conn):
    """Load genome_phenotypes.tsv as the growth_phenotypes_detailed table."""
    tsv_path = os.path.join(output_dir, "phenotypes", "genome_phenotypes.tsv")
    if not os.path.exists(tsv_path):
        print("No genome_phenotypes.tsv found, skipping growth_phenotypes_detailed table")
        return pd.DataFrame()

    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    df.to_sql("growth_phenotypes_detailed", conn, if_exists="replace", index=False)
    print(f"Built 'growth_phenotypes_detailed' table: {len(df)} rows")
    return df


# ---------------------------------------------------------------------------
# Table 9: gene_phenotypes
# ---------------------------------------------------------------------------
def build_gene_phenotypes_table(output_dir, conn):
    """Load gene_phenotypes.tsv as the gene_phenotypes table.

    Columns are loaded dynamically from whatever is in the TSV file.
    """
    tsv_path = os.path.join(output_dir, "phenotypes", "gene_phenotypes.tsv")
    if not os.path.exists(tsv_path):
        print("No gene_phenotypes.tsv found, skipping gene_phenotypes table")
        return pd.DataFrame()

    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    df.to_sql("gene_phenotypes", conn, if_exists="replace", index=False)
    print(f"Built 'gene_phenotypes' table: {len(df)} rows")
    return df


# ---------------------------------------------------------------------------
# Master function
# ---------------------------------------------------------------------------
def build_berdl_database(output_dir, db_path=None, reference_data_path="/data/"):
    """
    Build a complete BERDL SQLite database from pipeline output files.

    Args:
        output_dir: Path to pipeline output directory (e.g., .../Acinetobacter_baylyi_ADP1_RAST)
        db_path: Path for output SQLite database. Defaults to {output_dir}/berdl_tables.db
        reference_data_path: Path to reference data directory for ontology enrichment.
    """
    if db_path is None:
        db_path = os.path.join(output_dir, "berdl_tables.db")

    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    try:
        print(f"Building BERDL database: {db_path}")
        print(f"Input directory: {output_dir}")
        print("=" * 60)

        build_genome_table(output_dir, conn)
        build_genome_ani_table(output_dir, conn)
        build_genome_features_table(output_dir, conn)
        build_missing_functions_table(output_dir, conn)
        build_pan_genome_features_table(output_dir, conn)
        build_gene_reaction_data_table(output_dir, conn)
        build_genome_reactions_table(output_dir, conn)
        build_growth_phenotype_summary_table(output_dir, conn)
        build_growth_phenotypes_detailed_table(output_dir, conn)
        build_gene_phenotypes_table(output_dir, conn)

        conn.commit()
        print("=" * 60)
        print(f"Database complete: {db_path}")

        # Print summary
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        for (table_name,) in cursor.fetchall():
            cursor.execute(f"SELECT COUNT(*) FROM [{table_name}]")
            count = cursor.fetchone()[0]
            print(f"  {table_name}: {count} rows")
    finally:
        conn.close()

    # Generate ontology tables (operates on the closed database file)
    from KBDatalakeApps.KBDatalakeUtils import generate_ontology_tables
    generate_ontology_tables(
        input_database=db_path,
        reference_data_path=reference_data_path,
        source_tables=[
            "genome_features",
            "pan_genome_features",
            "genome_reactions",
            "growth_phenotypes_detailed",
        ],
    )


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python build_berdl_db.py <output_dir> [db_path]")
        sys.exit(1)
    output_dir = sys.argv[1]
    db_path = sys.argv[2] if len(sys.argv) > 2 else None
    build_berdl_database(output_dir, db_path)

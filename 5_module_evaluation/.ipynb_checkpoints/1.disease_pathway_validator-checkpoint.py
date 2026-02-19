#!/usr/bin/env python3
"""
Disease Pathway Pipeline - Modular Implementation

A generic, modular pipeline for disease-pathway enrichment analysis.

Workflow:
  0. Resolve disease term to stable IDs (DOID, MeSH, KEGG).
  1. Fetch associated genes using a selected gene provider (e.g., DisGeNET).
  2. Discover pathways, either by:
     a) 'enrichment': Performing enrichment analysis on the gene set (g:Profiler).
     b) 'direct': Fetching pre-curated pathways from a source (e.g., KEGG, Reactome).
  3. Compute enrichment of user-provided PPI modules against the discovered pathways.
  4. Rank modules by enrichment scores to find the most relevant ones.
"""

# python 1.disease_pathway_validator.py \
#   --disease "Herpes Simplex Virus infection" \
#   --modules /home/bbc8731/HSV/3_module_expansion/data/categories_methods \
#   --ppi /home/bbc8731/diseasemodulediscovery/tests/uniprot_ppi.csv \
#   --output-prefix HSV \
#   --pathway-source kegg reactome \
#   --n-methods-thresholds 1 2 3 4 5 6

import logging
logging.disable(logging.CRITICAL)

import argparse
import glob
import logging
import os
import sys
from functools import lru_cache
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
from gprofiler import GProfiler
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import requests
from pathlib import Path

# Import from local provider modules
from providers.disease_resolver import resolve
from providers.sources import GENE_SOURCES, PATHWAY_SOURCES

# --- Configuration ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Caching for I/O ---
# # I have added read_module_genes to remove cach for now
# def read_module_genes(module_file: str) -> Set[str]:
#     df = pd.read_csv(module_file, sep=None, engine='python', comment='#')

#     if 'n_methods' in df.columns:
#         print("_______Hihihihihi")
#         df = df.loc[df['n_methods'] > 5]

#     gene_columns = ['gene','Gene','GENE','gene_symbol','symbol','Symbol','name','node','protein']
#     for col in gene_columns:
#         if col in df.columns:
#             return set(df[col].dropna().astype(str))

#     if len(df.columns) == 2:
#         return set(df.iloc[:, 0].astype(str)) | set(df.iloc[:, 1].astype(str))

#     return set(df.iloc[:, 0].astype(str)) if not df.empty else set()

@lru_cache(maxsize=None)
def read_module_genes(module_file: str, n_methods_gt: int) -> Set[str]:
    """
    Load genes from a module file, with caching to avoid re-reading.
    This version correctly handles both node lists and edge lists.
    """
    try:
        df = pd.read_csv(module_file, sep=None, engine='python', comment='#') # Auto-detect separator, ignore comments
        if 'n_methods' in df.columns:
            df = df.loc[df['n_methods'] >= n_methods_gt, :]

        # 1. Look for a standard node/gene list column header.
        #    'name' is added to handle your second example robustly.
        gene_columns = ['gene', 'Gene', 'GENE', 'gene_symbol', 'symbol', 'Symbol', 'name', 'node', 'protein']
        for col in gene_columns:
            if col in df.columns:
                logger.debug(f"Reading gene list from column '{col}' in {module_file}")
                return set(df[col].dropna().astype(str))

        # 2. If no standard header is found, check for a two-column edge list format.
        #    This correctly handles your first example.
        if len(df.columns) == 2:
            logger.debug(f"Assuming two-column edge list format for {module_file}")
            col1_genes = set(df.iloc[:, 0].dropna().astype(str))
            col2_genes = set(df.iloc[:, 1].dropna().astype(str))
            # Combine genes from both columns
            all_genes = col1_genes.union(col2_genes)
            return all_genes

        # 3. As a last resort, fall back to using the first column.
        if not df.empty:
            logger.debug(f"Falling back to first column for {module_file}")
            return set(df.iloc[:, 0].dropna().astype(str))

    except Exception as e:
        logger.warning(f"Error loading module genes from {module_file}: {e}")
    
    return set()
    
# --- Core Pipeline Class ---
class DiseasePathwayPipeline:
    def __init__(self, output_prefix: str, base_dir: Path):
        self.output_prefix = output_prefix
        self.base_dir = base_dir
        self.universe_genes = None
        self.gp = GProfiler(return_dataframe=True, user_agent="DiseasePathwayPipeline/2.0")
    
    def discover_pathways_from_genes(self, gene_list: List[str]) -> pd.DataFrame:
        """(METHOD A) Find enriched pathways from a gene list using g:Profiler."""
        logger.info(f"Discovering pathways via gene set enrichment for {len(gene_list)} genes.")
        if not gene_list:
            logger.warning("Gene list is empty, cannot perform enrichment.")
            return pd.DataFrame()
        
        result = self.gp.profile(
            organism='hsapiens', query=gene_list,
            sources=['KEGG', 'REAC', 'WP'], user_threshold=0.05,
            significance_threshold_method='fdr', no_evidences=False
        )
        if result.empty: return pd.DataFrame()

        result['genes'] = result['intersections'].apply(lambda x: ';'.join(x))
        result = result.rename(columns={'native': 'pathway_id', 'source': 'database', 'name': 'name'})
        return result[['pathway_id', 'database', 'name', 'genes', 'p_value']]

    def compute_and_rank(
        self,
        pathways_df: pd.DataFrame,
        module_files: List[str],
        fdr_thr: float,
        rank_metric: str,
        ppi_universe: Set[str],
        n_methods_gt: int
    ):
        """Runs the enrichment and ranking steps of the pipeline."""
        if pathways_df.empty:
            logger.error("No pathways were discovered. Aborting enrichment and ranking.")
            return pd.DataFrame(), pd.DataFrame()

        logger.info(f"Starting enrichment against {len(pathways_df)} pathways...")
        enrichment_df = self._compute_enrichment(pathways_df, module_files, fdr_thr, ppi_universe, n_methods_gt)


        logger.info("Starting module ranking...")
        ranking_df = self._rank_modules(enrichment_df, rank_metric)
        
        return enrichment_df, ranking_df


    def _compute_enrichment(
        self,
        pathways_df: pd.DataFrame,
        module_files: List[str],
        fdr_thr: float,
        ppi_universe: Set[str],
        n_methods_gt: int
    ) -> pd.DataFrame:
    
        self.universe_genes = ppi_universe

        # self.universe_genes = self._build_dynamic_universe(pathways_df, module_files)
        results = []
        for module_file in module_files:
            p = Path(module_file)
            module_genes = read_module_genes(module_file, n_methods_gt)
            if not module_genes: continue

            try:
                category = p.parents[2].name
            except IndexError:
                category = "unknown_category"

            for _, p_row in pathways_df.iterrows():
                pathway_genes = set(p_row['genes'].split(';')) if pd.notna(p_row['genes']) else set()
                if not pathway_genes: continue
                
                enr = self._hypergeometric_test(module_genes, pathway_genes, self.universe_genes)
                results.append({
                    'module': p.relative_to(self.base_dir).as_posix(),
                    'pathway_id': p_row['pathway_id'],
                    'pathway_name': p_row['name'],
                    'pathway_database': p_row['database'],
                    'module_size': len(module_genes),
                    'pathway_size': len(pathway_genes),
                    'overlap': enr['overlap'],
                    'p_value': enr['p_value'],
                    'fold_enrichment': enr['fold_enrichment']
                })

        if not results: return pd.DataFrame()
        enrichment_df = pd.DataFrame(results)
        enrichment_df['fdr'] = multipletests(enrichment_df['p_value'], method='fdr_bh')[1]
        enrichment_df['significant'] = enrichment_df['fdr'] < fdr_thr
        
        # --- THIS IS THE NEW LINE ---
        # Sort the detailed results by FDR and then p-value for better readability.
        logger.info("Sorting enrichment results by significance (FDR, then p-value).")
        enrichment_df = enrichment_df.sort_values(by=['fdr', 'p_value'], ascending=True).reset_index(drop=True)
        # --- END OF NEW LINE ---
        
        out_file = f"data/{self.output_prefix}_module_enrichment.csv"
        enrichment_df.to_csv(out_file, index=False)
        logger.info(f"Enrichment results saved to {out_file}")
        return enrichment_df


    def _rank_modules(self, enrichment_df: pd.DataFrame, rank_metric: str) -> pd.DataFrame:
        if enrichment_df.empty: return pd.DataFrame()
        
        scores = []
        # --- NEW: Group by both category and module ---
        for (module), m_data in enrichment_df.groupby(['module']):
            log_p = -np.log10(m_data['p_value'].clip(lower=1e-300))
            
            weights = m_data['pathway_size']
            weighted_sumlog_score = (weights * log_p).sum() / weights.sum() if weights.sum() > 0 else 0
            mean_log_p_score = log_p.mean()
            max_enrichment_score = m_data['fold_enrichment'].max()
            
            if rank_metric == "weighted_sumlog":
                primary_score = weighted_sumlog_score
            elif rank_metric == "mean_log_p":
                primary_score = mean_log_p_score
            else:
                primary_score = max_enrichment_score

            scores.append({
                'module': module,
                'primary_ranking_score': primary_score,
                'rank_metric_used': rank_metric,
                'significant_pathways_count (FDR<0.05)': (m_data['fdr'] < 0.05).sum(),
                'score_weighted_sumlog': weighted_sumlog_score,
                'score_mean_log10_p': mean_log_p_score,
                'max_fold_enrichment': max_enrichment_score,
                'min_p_value': m_data['p_value'].min()
            })

        ranking_df = pd.DataFrame(scores)
        ranking_df = ranking_df.sort_values('primary_ranking_score', ascending=False).reset_index(drop=True)
        ranking_df['rank'] = ranking_df.index + 1
        
        # Reorder columns to put category first after rank
        cols_order = ['rank', 'module', 'primary_ranking_score', 'rank_metric_used', 
                      'significant_pathways_count (FDR<0.05)', 'score_weighted_sumlog', 
                      'score_mean_log10_p', 'max_fold_enrichment', 'min_p_value']
        ranking_df = ranking_df[cols_order]
        
        out_file = f"data/{self.output_prefix}_module_ranking.csv"
        ranking_df.to_csv(out_file, index=False)
        logger.info(f"Module ranking saved to {out_file}")
        return ranking_df
    
    def _build_dynamic_universe(self, ppi_universe):
        print("Universe size:", len(ppi_universe))
        return ppi_universe

    
    def _hypergeometric_test(self, mod_genes: Set, path_genes: Set, uni_genes: Set) -> Dict:
        mod_in_uni = mod_genes & uni_genes
        path_in_uni = path_genes & uni_genes
        overlap = len(mod_in_uni & path_in_uni)
        M, n, N = len(uni_genes), len(path_in_uni), len(mod_in_uni)
        if N == 0 or n == 0: return {'overlap': 0, 'p_value': 1.0, 'fold_enrichment': 0.0}
        p_val = hypergeom.sf(overlap - 1, M, n, N)
        expected = (N * n) / M
        fold = overlap / expected if expected > 0 else 0.0
        return {'overlap': overlap, 'p_value': p_val, 'fold_enrichment': fold}



def load_ppi_universe(ppi_file: str) -> Set[str]:
    df = pd.read_csv(ppi_file, sep=None, engine="python")
    if df.shape[1] >= 2:
        return set(df.iloc[:, 0].astype(str)) | set(df.iloc[:, 1].astype(str))
    else:
        return set(df.iloc[:, 0].astype(str))
    


# --- Main Execution ---
# (In disease_pathway_validator_modular.py)

def main():
    resolve.cache_clear()
    # read_module_genes.cache_clear()
    parser = argparse.ArgumentParser(
        description="Disease Pathway Enrichment Pipeline (Modular Version)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--disease", required=True, help="Disease name to analyze.")
    # parser.add_argument("--modules", nargs="+", required=True, help="PPI module files.")
    parser.add_argument("--modules", required=True, help="Path to results_full_PPI directory")

    parser.add_argument("--output-prefix", default="validation", help="Prefix for output files.")
    
    # --- MODIFICATION: Accept one or more pathway sources ---
    parser.add_argument(
        "--pathway-source", nargs='+', default=["enrichment"],
        choices=["enrichment"] + list(PATHWAY_SOURCES.keys()),
        help="One or more methods to discover pathways. 'enrichment' uses g:Profiler. Others use direct lookups."
    )
    # --- END MODIFICATION ---
    
    parser.add_argument("--gene-source", default="disgenet", choices=list(GENE_SOURCES.keys()), help="Database for gene-disease associations.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05, help="FDR threshold for significance.")
    parser.add_argument("--rank-metric", default="weighted_sumlog", choices=["weighted_sumlog", "mean_log_p", "max_enrichment"], help="Metric for ranking modules.")
    # Add back the disgenet threshold argument if you use the enrichment method
    parser.add_argument("--disgenet-threshold", type=float, default=0.3, help="Score threshold for DisGeNET.")
    parser.add_argument("--ppi", required=True, help="Path to PPI file (edge list or node list) used for module detection")
    parser.add_argument(
        "--n-methods-thresholds",
        nargs="+",
        type=int,
        default=[1, 2, 3, 4, 5, 6],
        help="Run one analysis per cutoff, filtering module genes with n_methods > cutoff."
    )

    args = parser.parse_args()
    ppi_universe = load_ppi_universe(args.ppi)
    print("PPI universe size:", len(ppi_universe))

    base_dir = Path(args.modules)

    # --- Pipeline Execution ---
    resolved_ids = resolve(args.disease)
    if not resolved_ids: sys.exit("Could not resolve disease term. Exiting.")

    pipeline = DiseasePathwayPipeline(output_prefix=args.output_prefix, base_dir=base_dir)
    all_pathway_dfs = []

    # Handle the 'enrichment' case if requested
    if "enrichment" in args.pathway_source:
        logger.info(f"Using ENRICHMENT method with gene source: {args.gene_source}")
        gene_provider = GENE_SOURCES[args.gene_source]
        disease_genes = gene_provider.get_genes(resolved_ids, score_threshold=args.disgenet_threshold)
        pathways_df = pipeline.discover_pathways_from_genes(disease_genes)
        if not pathways_df.empty: all_pathway_dfs.append(pathways_df)

    # Handle all requested direct provider cases
    direct_sources = [source for source in args.pathway_source if source != "enrichment"]
    for source_name in direct_sources:
        logger.info(f"Using DIRECT pathway source: {source_name}")
        if source_name in PATHWAY_SOURCES:
            pathway_provider = PATHWAY_SOURCES[source_name]
            # --- THIS IS THE MODIFIED LINE ---
            # Pass the original disease name to the provider
            pathways_df = pathway_provider.get_pathways(resolved_ids, disease_name=args.disease)
            # --- END OF MODIFICATION ---
            if not pathways_df.empty:
                all_pathway_dfs.append(pathways_df)
        else:
            logger.warning(f"Unknown pathway source '{source_name}' requested. Skipping.")
    
    final_pathways_df = pd.concat(all_pathway_dfs, ignore_index=True).drop_duplicates(subset=['pathway_id']).reset_index(drop=True)

    out_file = f"data/{args.output_prefix}_disease_specific_pathways.csv"
    final_pathways_df.to_csv(out_file, index=False)
    logger.info(f"Discovered a combined total of {len(final_pathways_df)} unique pathways, saved to {out_file}")
    
    # module_files = args.modules
    module_files = list(base_dir.rglob("*/consensus/uniprot_ppi.tsv"))   # consenses of methods, only ppi

    print("Modules found:", len(module_files))
    print("Pathways found:", len(final_pathways_df))

    thresholds = sorted(set(args.n_methods_thresholds))
    print(f"n_methods cutoffs to run (n_methods > cutoff): {thresholds}")

    for cutoff in thresholds:
        run_prefix = f"{args.output_prefix}_n_methods_gt_{cutoff}"
        run_pipeline = DiseasePathwayPipeline(output_prefix=run_prefix, base_dir=base_dir)
        print(f"\nRunning with filter: n_methods > {cutoff}")
        enrichment_df, ranking_df = run_pipeline.compute_and_rank(
            final_pathways_df,
            module_files,
            args.fdr_threshold,
            args.rank_metric,
            ppi_universe,
            cutoff
        )
        if not ranking_df.empty:
            top_module = ranking_df.iloc[0]
            print(f"Top module (n_methods > {cutoff}): {top_module['module']}")
        else:
            print(f"No ranked modules for n_methods > {cutoff}")

    # --- Final Summary ---
    print("\n" + "="*60 + "\nPIPELINE SUMMARY\n" + "="*60)
    print(f"Disease: '{args.disease}' (Resolved to: {resolved_ids})")
    print(f"Pathway Sources Used: {', '.join(s.upper() for s in args.pathway_source)}")
    print(f"Total Unique Pathways Discovered: {len(final_pathways_df)}")
    print(f"Modules Analyzed: {len(module_files)}")
    print("="*60)

if __name__ == "__main__":
    main()

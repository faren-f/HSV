# (In providers/sources.py)

import logging
from typing import Protocol, List, Dict
import pandas as pd
import requests
import re
import time

logger = logging.getLogger(__name__)
USER_AGENT = "DiseasePathwayPipeline/2.0 (Modular)"
REQUEST_TIMEOUT = 30

# --- Interfaces (Protocols) ---
# (No changes needed here)
class GeneSource(Protocol):
    """Interface for providers that return a list of disease-associated genes."""
    def get_genes(self, ids: Dict[str, str], **kwargs) -> List[str]: ...

class PathwaySource(Protocol):
    """Interface for providers that return a DataFrame of pre-defined disease pathways."""
    def get_pathways(self, ids: Dict[str, str], **kwargs) -> pd.DataFrame: ...

# --- Helper Functions ---

# --- THIS IS THE CORRECTED HELPER FUNCTION ---
def _api_get(url, params=None, extra_headers=None):
    """
    Generic GET request helper with error handling and flexible headers.
    """
    headers = {"User-Agent": USER_AGENT}
    if extra_headers:
        headers.update(extra_headers) # Merge in any additional headers

    try:
        response = requests.get(url, params=params, timeout=REQUEST_TIMEOUT, headers=headers)
        response.raise_for_status()
        return response
    except requests.RequestException as e:
        logger.warning(f"API request failed for {url}: {e}")
        return None
# --- END OF CORRECTED HELPER FUNCTION ---

# --- Concrete Gene Source Implementations ---
# (No changes needed here)
class DisGeNetGeneSource(GeneSource):
    """Fetches disease-associated genes from DisGeNET using MeSH or DOID."""
    def get_genes(self, ids: Dict[str, str], score_threshold: float = 0.3, max_genes: int = 500) -> List[str]:
        disease_id = ids.get("MESH") or ids.get("DOID")
        if not disease_id:
            logger.warning("DisGeNET provider requires a MESH or DOID, none found.")
            return []
        logger.info(f"Querying DisGeNET with ID: {disease_id}")
        gda_url = "https://www.disgenet.org/api/gda/disease"
        params = {"disease": disease_id, "format": "json", "limit": max_genes}
        response = _api_get(gda_url, params=params)
        if not response:
            return []

        print("Status code:", response.status_code)
        print("Response text:", response.text[:500])

        associations = response.json()
        genes = [
            assoc['gene_symbol'] for assoc in associations
            if float(assoc.get('score', 0)) >= score_threshold and assoc.get('gene_symbol')
        ]
        logger.info(f"DisGeNET found {len(genes)} genes with score >= {score_threshold}")
        return genes

# --- Concrete Pathway Source Implementations ---

# (No changes needed in KeggPathwaySource)
# (Helper functions for KeggPathwaySource)
def batch_list(input_list, batch_size):
    """Yield successive n-sized chunks from a list."""
    for i in range(0, len(input_list), batch_size):
        yield input_list[i:i + batch_size]

class KeggPathwaySource(PathwaySource):
    """
    [PERFORMANCE-TUNED VERSION]
    Fetches disease pathways from KEGG using batched calls and a
    precise regular expression to extract ONLY valid UniProt IDs.
    """
    def get_pathways(self, ids: Dict[str, str], **kwargs) -> pd.DataFrame:
        print(">>> KEGG PROVIDER CALLED")   

        kegg_id = ids.get("KEGG")
        if not kegg_id:
            logger.warning("KEGG pathway source requires a KEGG Disease ID, none found.")
            return pd.DataFrame()
        
        logger.info(f"Querying KEGG for pathways associated with disease {kegg_id}")
        link_res = _api_get(f"http://rest.kegg.jp/link/pathway/{kegg_id}")
        if not link_res or not link_res.text.strip():
            logger.warning(f"No KEGG pathways found linked to {kegg_id}")
            return pd.DataFrame()

        pathway_ids_with_prefix = {line.split('\t')[1] for line in link_res.text.strip().split('\n') if 'path:hsa' in line}
        
        if not pathway_ids_with_prefix:
            logger.warning(f"Found related KEGG entries for {kegg_id}, but none were human pathways.")
            return pd.DataFrame()

        logger.info(f"Found {len(pathway_ids_with_prefix)} human pathway IDs: {', '.join(pathway_ids_with_prefix)}")
        
        pathways = []
        for pid_with_prefix in pathway_ids_with_prefix:
            pw_res = _api_get(f"http://rest.kegg.jp/get/{pid_with_prefix}")
            if not pw_res: continue
            
            name = [l.replace("NAME", "").strip() for l in pw_res.text.split('\n') if l.startswith("NAME")][0]
            pathway_id = pid_with_prefix.replace('path:', '')
            gene_res = _api_get(f"http://rest.kegg.jp/link/hsa/{pathway_id}")
            
            if gene_res and gene_res.text.strip():
                kegg_gene_ids = [line.split('\t')[1] for line in gene_res.text.strip().split('\n')]
                
                uniprot_ids = set()
                logger.info(f"Batch-fetching UniProt IDs from KEGG for {len(kegg_gene_ids)} genes in {pathway_id}...")
                
                for gene_batch in batch_list(kegg_gene_ids, 10):
                    batch_query = "+".join(gene_batch)
                    batch_entry_res = _api_get(f"http://rest.kegg.jp/get/{batch_query}")
                    if batch_entry_res:
                        for entry_text in batch_entry_res.text.strip().split('///\n'):
                            if not entry_text: continue
                            uniprot_line_match = re.search(r"UniProt:\s+(.*)", entry_text)
                            if uniprot_line_match:
                                uniprot_line = uniprot_line_match.group(1)
                                found_ids = re.findall(r"\b[A-Z0-9]{6,10}\b", uniprot_line)
                                uniprot_ids.update(found_ids)
                    time.sleep(0.1)
                
                if uniprot_ids:
                    logger.info(f"Found {len(uniprot_ids)} valid UniProt IDs for pathway {pathway_id}.")
                    pathways.append({
                        'pathway_id': pathway_id, 'database': 'KEGG', 
                        'name': name, 'genes': ';'.join(sorted(list(uniprot_ids)))
                    })
        
        logger.info(f"KEGG source successfully retrieved {len(pathways)} direct pathways with gene sets.")
        return pd.DataFrame(pathways)

class ReactomePathwaySource(PathwaySource):
    """
    [IMPROVED VERSION]
    Fetches disease pathways from Reactome using a text-based search
    to find relevant pathways even without a direct DOID mapping.
    """
    def get_pathways(self, ids: Dict[str, str], **kwargs) -> pd.DataFrame:
        # Use the original disease term from the main script for the search.
        # This requires a small change in the main script to pass the term.
        # For now, we'll fall back to using what we can from the IDs.
        print(">>> REACTOME PROVIDER CALLED")   

        disease_name = kwargs.get("disease_name", list(ids.keys())[0] if ids else "disease")

        logger.info(f"Querying Reactome with search term: '{disease_name}'")

        search_url = "https://reactome.org/ContentService/search/query"
        params = {'query': disease_name, 'species': 'Homo sapiens', 'types': 'Pathway', 'format': 'json'}
        
        response = _api_get(search_url, params=params)
        
        if not response or not response.json().get('results'):
            logger.warning(f"Reactome text search found no pathways for '{disease_name}'")
            return pd.DataFrame()
            
        search_results = response.json()['results']
        
        # --- Intelligent Filtering of Search Results ---
        # Generate keywords to ensure we only get relevant pathways.
        disease_lower = disease_name.lower()
        keywords = {disease_lower}
        if 'herpes' in disease_lower or 'hsv' in disease_lower:
            keywords.update(['herpes', 'hsv', 'viral infection'])
        if 'alzheimer' in disease_lower:
            keywords.update(['alzheimer', 'neurodegeneration'])
        
        relevant_pathways = []
        for pathway_data in search_results:
            pathway_name_lower = pathway_data.get('name', '').lower()
            # Check if any of our keywords are in the found pathway's name
            if any(key in pathway_name_lower for key in keywords):
                relevant_pathways.append(pathway_data)

        if not relevant_pathways:
            logger.warning(f"Reactome search returned results, but none were deemed relevant to '{disease_name}' after filtering.")
            return pd.DataFrame()

        logger.info(f"Found {len(relevant_pathways)} potentially relevant pathways in Reactome. Fetching genes...")

        pathways = []
        for pw_data in relevant_pathways:
            stid = pw_data['stId']
            # Get genes for each pathway
            proj_url = f"https://reactome.org/ContentService/data/participants/{stid}/accessions"
            proj_res = _api_get(proj_url, extra_headers={"accept": "text/plain"})
            
            if proj_res and proj_res.text.strip():
                uniprot_ids = set(proj_res.text.strip().split('\n'))
                if uniprot_ids:
                    pathways.append({
                        'pathway_id': stid, 'database': 'REAC',
                        'name': pw_data['name'], 'genes': ';'.join(sorted(list(uniprot_ids)))
                    })
        
        logger.info(f"Reactome source successfully retrieved {len(pathways)} direct pathways with gene sets.")
        return pd.DataFrame(pathways)

# --- Provider Registries ---
# (No changes needed here)
GENE_SOURCES: Dict[str, GeneSource] = {
    "disgenet": DisGeNetGeneSource(),
}

PATHWAY_SOURCES: Dict[str, PathwaySource] = {
    "kegg": KeggPathwaySource(),
    "reactome": ReactomePathwaySource(),
}
# providers/disease_resolver.py

import logging
import re
from functools import lru_cache

import requests

logger = logging.getLogger(__name__)
USER_AGENT = "DiseasePathwayPipeline/2.0 (Modular)"
REQUEST_TIMEOUT = 15

@lru_cache(maxsize=32)
def resolve(term: str) -> dict[str, str]:
    """
    Map an arbitrary free-text disease term to a dictionary of stable IDs.

    This function queries multiple sources to find the best possible match and
    returns a dictionary containing available identifiers like DOID, MeSH, and KEGG.

    Args:
        term: The free-text disease term (e.g., "herpes simplex virus infection").

    Returns:
        A dictionary of found IDs, e.g.,
        {"DOID": "DOID:10563", "MESH": "D006561", "KEGG": "H00371"}
    """
    ids = {}
    logger.info(f"Resolving term '{term}' to stable disease IDs...")

    # 1) EBI OLS (Disease Ontology) -> DOID, MeSH
    try:
        ols_url = "https://www.ebi.ac.uk/ols/api/search"
        params = {
            "q": term,
            "ontology": "doid",
            "exact": "true",
            "fieldList": "obo_id,label,annotation",
            "rows": 1
        }
        response = requests.get(ols_url, params=params, timeout=REQUEST_TIMEOUT, headers={"User-Agent": USER_AGENT})
        response.raise_for_status()
        ols_data = response.json()

        if ols_data["response"]["docs"]:
            doc = ols_data["response"]["docs"][0]
            ids["DOID"] = doc.get("obo_id", "")
            if "annotation" in doc and "hasDbXref" in doc["annotation"]:
                # Find the first MeSH ID in the cross-references
                mesh_id = next((x for x in doc["annotation"]["hasDbXref"] if x.startswith("MESH:")), "")
                ids["MESH"] = mesh_id.replace("MESH:", "")
            logger.info(f"  > Found via EBI OLS: DOID={ids.get('DOID')}, MESH={ids.get('MESH')}")
    except requests.RequestException as e:
        logger.warning(f"Could not query EBI OLS API: {e}")

    # 2) KEGG DISEASE (regex match) -> KEGG
    try:
        kegg_url = f"http://rest.kegg.jp/find/disease/{term}"
        response = requests.get(kegg_url, timeout=REQUEST_TIMEOUT, headers={"User-Agent": USER_AGENT})
        response.raise_for_status()
        # Find the first human disease ID (H-number)
        match = re.search(r"(H\d+)", response.text)
        if match:
            ids["KEGG"] = match.group(1)
            logger.info(f"  > Found via KEGG: KEGG={ids.get('KEGG')}")
    except requests.RequestException as e:
        logger.warning(f"Could not query KEGG API: {e}")

    if not ids:
        logger.error(f"Could not resolve '{term}' to any known disease ID.")

    return {k: v for k, v in ids.items() if v} # Return only non-empty IDs
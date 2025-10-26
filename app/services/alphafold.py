"""
AlphaFold Service

Fetch AlphaFold predictions and download CIF files
"""
import requests
from pathlib import Path
from typing import List, Dict, Optional
from app.util import backoff


# Project root: protein-service/
BASE_DIR = Path(__file__).resolve().parents[2]
CACHE = (BASE_DIR / "data" / "alphafold")
CACHE.mkdir(parents=True, exist_ok=True)


SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "protein-service/1.0"})


@backoff()
def _get(url, timeout=30):
    """HTTP GET with retry logic"""
    return SESSION.get(url, timeout=timeout)


def fetch_alphafold_models(uniprot_id: str) -> List[Dict]:
    """
    Fetch AlphaFold model metadata from API
    
    Args:
        uniprot_id: UniProt accession ID
        
    Returns:
        List of model metadata dicts
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        r = _get(url)
    except requests.exceptions.Timeout:
        return []
    except requests.exceptions.RequestException:
        return []
    
    if r.status_code != 200:
        return []
    
    try:
        arr = r.json()
    except Exception:
        return []
    
    if not isinstance(arr, list):
        return []
    
    out = []
    for m in arr:
        out.append({
            "uniprot_id": m.get("uniprotAccession") or m.get("uniprot_id"),
            "pdbUrl": m.get("pdbUrl"),
            "cifUrl": m.get("cifUrl"),
            "pae_image_url": m.get("paeImageUrl"),
            "plddt": m.get("plddt") or m.get("pLDDT"),
        })
    return out


def download_alphafold_cif(uniprot_id: str, cif_url: str) -> Optional[str]:
    """
    Download AlphaFold CIF file
    
    Args:
        uniprot_id: UniProt accession ID (for filename)
        cif_url: URL to CIF file
        
    Returns:
        Path to downloaded file or None if failed
    """
    fn = CACHE / f"AF-{uniprot_id}-F1-model_v4.cif"
    
    # Use cached file if exists
    if fn.exists() and fn.stat().st_size > 0:
        return str(fn)
    
    # Download
    try:
        r = _get(cif_url)
    except Exception:
        return None
    
    if r.status_code == 200 and r.content:
        fn.write_bytes(r.content)
        return str(fn)
    
    return None


def fetch_alphafold_with_cif(uniprot_id: str) -> tuple[List[Dict], List[str]]:
    """
    Fetch AlphaFold models and download CIF files
    
    Args:
        uniprot_id: UniProt accession ID
        
    Returns:
        Tuple of (model_metadata_list, cif_path_list)
    """
    models = fetch_alphafold_models(uniprot_id)
    cif_paths = []
    
    for m in models:
        cif_url = m.get("cifUrl")
        if cif_url:
            path = download_alphafold_cif(uniprot_id, cif_url)
            if path:
                cif_paths.append(path)
    
    return models, cif_paths

"""
PDB (RCSB) Service

Download mmCIF structure files from RCSB PDB
"""
import requests
from pathlib import Path
from typing import List, Optional


# Project root: protein-service/
BASE_DIR = Path(__file__).resolve().parents[2]  # app/services/pdb.py -> app/ -> protein-service/
CACHE = (BASE_DIR / "data" / "pdb")
CACHE.mkdir(parents=True, exist_ok=True)


def _get(url, timeout=20):
    """HTTP GET with timeout"""
    return requests.get(url, timeout=timeout)


def download_mmCIF(pdb_id: str) -> Optional[str]:
    """
    Download mmCIF file for a PDB ID
    
    Args:
        pdb_id: PDB ID (4 characters)
        
    Returns:
        Path to downloaded file or None if failed
    """
    fn = CACHE / f"{pdb_id}.cif"
    
    # Use cached file if exists
    if fn.exists() and fn.stat().st_size > 0:
        return str(fn)
    
    # Download from RCSB
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        r = _get(url)
    except Exception:
        return None
    
    if r.status_code == 200 and r.content:
        fn.write_bytes(r.content)
        return str(fn)
    
    return None


def download_mmCIF_batch(pdb_ids: List[str]) -> List[str]:
    """
    Download multiple mmCIF files
    
    Args:
        pdb_ids: List of PDB IDs
        
    Returns:
        List of successfully downloaded file paths
    """
    out = []
    for pid in pdb_ids:
        p = download_mmCIF(pid)
        if p:
            out.append(p)
    return out

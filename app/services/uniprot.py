"""
UniProt API Service

Fetch protein information and PDB cross-references from UniProt REST API
"""
import requests
from typing import Dict, List, Optional, Tuple
from app.util import backoff


SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "protein-service/1.0"})


@backoff()
def _get(url, timeout=15):
    """HTTP GET with retry logic"""
    return SESSION.get(url, timeout=timeout)


def fetch_uniprot_core(uniprot_id: str) -> Optional[Dict]:
    """
    Fetch core protein information from UniProt
    
    Args:
        uniprot_id: UniProt accession ID
        
    Returns:
        Dict with id, name, length, organism or None if not found
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        r = _get(url)
    except Exception:
        return None
        
    if r.status_code != 200:
        return None
        
    try:
        j = r.json()
    except Exception:
        return None
    
    # Extract protein name
    name = None
    try:
        name = j.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
    except Exception:
        pass
    
    # Extract sequence length
    length = j.get("sequence", {}).get("length")
    
    # Extract organism
    organism = j.get("organism", {}).get("scientificName")
    
    return {
        "id": uniprot_id,
        "name": name,
        "length": length,
        "organism": organism
    }


def fetch_pdb_ids(
    uniprot_id: str,
    method: str = "",
    max_pdb: int = 8,
    max_resolution: float | None = 3.5,
) -> List[str]:
    """
    Fetch PDB IDs from UniProt cross-references
    
    Args:
        uniprot_id: UniProt accession ID
        method: Experimental method (X-ray, EM, NMR) - case insensitive
        max_pdb: Maximum number of PDB IDs to return
        max_resolution: Maximum resolution in Angstroms (None for no limit)
        
    Returns:
        List of PDB IDs sorted by resolution (best first)
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        r = _get(url)
    except Exception:
        return []
        
    if r.status_code != 200:
        return []
        
    try:
        j = r.json()
    except Exception:
        return []
    
    xrefs = j.get("uniProtKBCrossReferences", [])
    
    items: List[Tuple[str, float]] = []
    
    for x in xrefs:
        if x.get("database") != "PDB":
            continue
            
        pid = x.get("id")
        if not pid:
            continue
        
        # Extract properties
        props = {p["key"]: p.get("value") for p in x.get("properties", [])}
        
        # Check method
        meth = (props.get("Method") or props.get("method") or "").lower()
        if method and meth and meth != method.lower():
            continue
        
        # Extract resolution
        res_str = props.get("Resolution") or props.get("resolution") or ""
        try:
            res = float(res_str.replace("Å", "").replace("A", "").strip())
        except Exception:
            res = 99.9
        
        # Filter by resolution
        if (max_resolution is not None) and (res > max_resolution):
            continue
        
        items.append((pid[:4].upper(), res))
    
    # Sort by resolution and deduplicate
    seen = set()
    unique = []
    for pid, res in sorted(items, key=lambda t: t[1]):
        if len(pid) == 4 and pid.isalnum() and pid not in seen:
            seen.add(pid)
            unique.append(pid)
        if len(unique) >= max_pdb:
            break
    
    return unique

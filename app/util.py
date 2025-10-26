"""
Utility functions for retry logic and deduplication
"""
import time
import random
import functools


def backoff(retries=3, base=0.5, factor=2.0):
    """
    Retry decorator with exponential backoff
    
    Args:
        retries: Number of retry attempts
        base: Initial delay in seconds
        factor: Multiplier for each retry
    """
    def deco(fn):
        @functools.wraps(fn)
        def wrap(*a, **kw):
            delay = base
            for i in range(retries):
                try:
                    return fn(*a, **kw)
                except Exception as e:
                    if i == retries - 1:
                        raise
                    time.sleep(delay + random.uniform(0, 0.25))
                    delay *= factor
        return wrap
    return deco


def uniq(seq):
    """
    Remove duplicates while preserving order
    
    Args:
        seq: Input sequence
        
    Returns:
        List with duplicates removed
    """
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

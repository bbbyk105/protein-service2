#!/usr/bin/env python3
"""
FastAPI テストスクリプト

研究室のサンプル（P17456, C6H0Y9）で動作確認
"""
import requests
import json
import time


API_URL = "http://localhost:8000"


def test_health():
    """ヘルスチェック"""
    print("=" * 60)
    print("TEST 1: Health Check")
    print("=" * 60)
    
    response = requests.get(f"{API_URL}/health")
    print(f"Status: {response.status_code}")
    print(f"Response: {json.dumps(response.json(), indent=2)}")
    print()
    
    return response.status_code == 200


def test_analyze_single():
    """単一ID解析"""
    print("=" * 60)
    print("TEST 2: Single UniProt ID Analysis")
    print("=" * 60)
    
    payload = {
        "uniprot_ids": "P17456",
        "method": "X-ray",
        "seq_ratio": 20
    }
    
    print(f"Request: {json.dumps(payload, indent=2)}")
    print()
    
    start = time.time()
    response = requests.post(f"{API_URL}/analyze", json=payload)
    elapsed = time.time() - start
    
    print(f"Status: {response.status_code}")
    print(f"Time: {elapsed:.2f}s")
    
    if response.status_code == 200:
        data = response.json()
        print(f"\nResults:")
        print(f"  - Entries: {len(data['results'])}")
        print(f"  - Warnings: {len(data['warnings'])}")
        
        for r in data['results']:
            print(f"\n  UniProt: {r['uniprot_id']}")
            print(f"  Status: {r['status']}")
            if 'core' in r:
                print(f"  Name: {r['core'].get('name')}")
                print(f"  Organism: {r['core'].get('organism')}")
            if 'analysis' in r:
                kpi = r['analysis'].get('kpi', {})
                print(f"  cis count: {kpi.get('cis_count')}")
                print(f"  midrange fraction: {kpi.get('midrange_dist_fraction')}")
        
        if data['warnings']:
            print(f"\n  Warnings:")
            for w in data['warnings']:
                print(f"    - {w}")
    else:
        print(f"Error: {response.text}")
    
    print()
    return response.status_code == 200


def test_analyze_multiple():
    """複数ID解析"""
    print("=" * 60)
    print("TEST 3: Multiple UniProt IDs Analysis")
    print("=" * 60)
    
    payload = {
        "uniprot_ids": ["P17456", "C6H0Y9"],
        "method": "X-ray",
        "seq_ratio": 20
    }
    
    print(f"Request: {json.dumps(payload, indent=2)}")
    print()
    
    start = time.time()
    response = requests.post(f"{API_URL}/analyze", json=payload)
    elapsed = time.time() - start
    
    print(f"Status: {response.status_code}")
    print(f"Time: {elapsed:.2f}s")
    
    if response.status_code == 200:
        data = response.json()
        print(f"\nResults:")
        print(f"  - Entries: {len(data['results'])}")
        print(f"  - Warnings: {len(data['warnings'])}")
        
        for r in data['results']:
            print(f"\n  [{r['uniprot_id']}]")
            print(f"    Status: {r['status']}")
            if 'analysis' in r and 'kpi' in r['analysis']:
                kpi = r['analysis']['kpi']
                print(f"    cis: {kpi.get('cis_count')}")
                print(f"    structures: {kpi.get('total_structures')}")
        
        if data['warnings']:
            print(f"\n  Warnings:")
            for w in data['warnings']:
                print(f"    - {w}")
    else:
        print(f"Error: {response.text}")
    
    print()
    return response.status_code == 200


def main():
    print("\n🧪 FastAPI Test Suite")
    print("=" * 60)
    print()
    
    try:
        # Test 1: Health
        if not test_health():
            print("❌ Health check failed")
            return
        
        # Test 2: Single ID
        if not test_analyze_single():
            print("❌ Single ID analysis failed")
            return
        
        # Test 3: Multiple IDs
        if not test_analyze_multiple():
            print("❌ Multiple IDs analysis failed")
            return
        
        print("=" * 60)
        print("✅ All tests passed!")
        print("=" * 60)
        
    except requests.exceptions.ConnectionError:
        print("❌ Connection Error: FastAPIが起動していません")
        print("   以下を実行してください:")
        print("   ./start.sh")
    except Exception as e:
        print(f"❌ Error: {e}")


if __name__ == "__main__":
    main()

import requests
import json

# Base URL for the Preprocessing Controller
BASE_URL = "http://127.0.0.1:8888/api/v1/preprocessing"

TEST_PROJECT_ID = "p_1767495922_c912e"
TEST_DATASET_ID = "ds_20260103_220522_0f5c3"



def test_qc_calculate():
    """
    Step 1: Test QC Calculation (GET/Compute metrics)
    """
    print("1. QC Calculation")

    url = f"{BASE_URL}/qc/calculate"
    payload = {
        "project_id": TEST_PROJECT_ID,
        "dataset_id": TEST_DATASET_ID,
        "organism": "Human"
    }

    print(f"POST {url}")
    try:
        response = requests.post(url, json=payload)

        if response.status_code == 200:
            print("[Success] QC Metrics Calculated.")
            print(json.dumps(response.json(), indent=2))
        else:
            print(f"[Failed] Status: {response.status_code}")
            print(response.text)

    except requests.exceptions.ConnectionError:
        print("[Error] Cannot connect to server.")


def test_qc_filter():
    """
    Step 2: Test QC Filtering (Apply Filter -> Create Snapshot)
    """
    print("2. QC Filter (Create Snapshot)")

    url = f"{BASE_URL}/qc/filter"

    # Payload matching FilterQCRequest DTO
    payload = {
        "project_id": TEST_PROJECT_ID,
        "dataset_id": TEST_DATASET_ID,
        "min_genes": 200,  # Filter cells with < 200 genes
        "min_cells": 3,  # Filter genes present in < 3 cells
        "max_counts": 50000,  # Optional: Max counts
        "pct_mt_max": 20.0,  # Optional: Max 20% Mitochondrial content
        "pct_hb_max": 5.0  # Optional: Max 5% Hemoglobin
    }

    print(f"POST {url}")
    print("Payload Config:", json.dumps(payload, indent=2))

    try:
        response = requests.post(url, json=payload)

        # We expect 201 Created
        if response.status_code == 201:
            data = response.json()
            print("\n[Success] Filtering Applied & Snapshot Created!")
            print(f"Snapshot ID:       {data['snapshot_id']}")
            print(f"Saved Path:        {data['snapshot_path']}")
            print(f"Cells Remaining:   {data['n_obs_remaining']}")
            print(f"Genes Remaining:   {data['n_vars_remaining']}")

            # Validation check
            if data['n_obs_remaining'] == 0:
                print("[Warning] 0 cells remaining! Filters might be too strict.")
        else:
            print(f"\n[Failed] Status Code: {response.status_code}")
            print("Error Details:", response.text)

    except requests.exceptions.ConnectionError:
        print("[Error] Cannot connect to server. Is uvicorn running?")


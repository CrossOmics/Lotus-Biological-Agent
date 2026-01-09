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
        "min_genes": 1000,  # Filter cells with < 200 genes
        "min_cells": 3,  # Filter genes present in < 3 cells
        "max_counts": 2500,  # Optional: Max counts
        "pct_mt_max": 2.0,  # Optional: Max 20% Mitochondrial content
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


def test_hvg_process():
    """
    Step 3: Test HVG Selection & Scaling.

    This test verifies:
    1. Automatic discovery of the latest 'QC Filtered' snapshot (since we don't provide source_snapshot_id).
    2. Correct execution of Normalize -> Log1p -> HVG -> Scale.
    3. Generation of the dispersion plot and a new snapshot.
    """
    print("HVG Selection & Scaling")

    url = f"{BASE_URL}/hvg"

    # Payload matching RunHVGRequest DTO
    # Note: 'source_snapshot_id' is deliberately OMITTED to test auto-discovery logic.
    payload = {
        "project_id": TEST_PROJECT_ID,
        "dataset_id": TEST_DATASET_ID,
        "n_top_genes": 2000,
        "flavor": "seurat",
        "target_sum": 1e4
    }

    print(f"POST {url}")
    print("Payload Config (Auto-find source snapshot):")
    print(json.dumps(payload, indent=2))

    try:
        response = requests.post(url, json=payload)

        # We expect 201 Created
        if response.status_code == 201:
            data = response.json()
            print("\n[Success] HVG & Scaling Completed!")
            print("-" * 30)
            print(f"Snapshot ID:       {data['snapshot_id']}")
            print(f"Saved Path:        {data['snapshot_path']}")
            print(f"HVG Plot Path:     {data['hvg_plot_path']}")
            print(f"Genes Found:       {data['n_genes_found']}")
            print(f"Message:           {data['msg']}")

            # Logic Check
            if data['n_genes_found'] != payload['n_top_genes']:
                print(f"\n[Warning] Requested {payload['n_top_genes']} genes but got {data['n_genes_found']}.")

        elif response.status_code == 404:
            print(f"\n[Failed] 404 Not Found")
            print("Server Message:", response.json().get('detail'))
            print("Tip: Did you run 'test_qc_filter' first to generate a QC snapshot?")

        elif response.status_code == 500:
            print(f"\n[Failed] 500 Internal Error")
            print("Error Details:", response.text)

        else:
            print(f"\n[Failed] Status Code: {response.status_code}")
            print("Response:", response.text)

    except requests.exceptions.ConnectionError:
        print(f"\n[Error] Could not connect to backend at {url}")
        print("Is the server running on port 8888?")

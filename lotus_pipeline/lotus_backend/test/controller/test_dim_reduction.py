import requests
import json

# Base URL for the Dimensionality Reduction Controller
BASE_URL = "http://127.0.0.1:8888/api/v1/dim-reduction"

TEST_PROJECT_ID = "p_1767495922_c912e"
TEST_DATASET_ID = "ds_20260103_220522_0f5c3"


def test_pca_analysis():
    """
    Step 1: Test PCA Analysis
    Performs Principal Component Analysis and generates Elbow Plot data.
    """
    print("1. PCA Analysis")
    url = f"{BASE_URL}/pca"

    # Payload matching RunPCARequest DTO
    payload = {
        "project_id": TEST_PROJECT_ID,
        "dataset_id": TEST_DATASET_ID,
        "n_comps": 50,  # Number of principal components to compute
        "svd_solver": "arpack"  # Optional: PCA solver method
    }

    print(f"POST {url}")
    print("Payload Config:", json.dumps(payload, indent=2))

    try:
        response = requests.post(url, json=payload)

        # We expect 201 Created
        if response.status_code == 201:
            data = response.json()
            print("\n[Success] PCA Analysis Completed!")
            print("-" * 30)
            print(f"Snapshot ID: {data['snapshot_id']}")
            print(f"Saved Path: {data['snapshot_path']}")
            print(f"pca_scatter_path: {data['pca_scatter_path']}")
            print(f"variance_plot_data: {data['variance_plot_data']}")

        elif response.status_code == 404:
            print(f"\n[Failed] 404 Not Found")
            print("Server Message:", response.json().get('detail'))
            print("Tip: Did you run HVG selection first to generate the required snapshot?")

        elif response.status_code == 500:
            print(f"\n[Failed] 500 Internal Error")
            print("Error Details:", response.text)

        else:
            print(f"\n[Failed] Status Code: {response.status_code}")
            print("Response:", response.text)

    except requests.exceptions.ConnectionError:
        print(f"\n[Error] Could not connect to backend at {url}")
        print("Is the server running on port 8888?")


def test_build_neighbors():
    """
    Step 2: Test Neighborhood Graph Construction
    Builds k-nearest neighbor graph based on PCA results.
    User should select n_pcs based on the Elbow Plot from PCA analysis.
    """
    print("\n2. Build Neighborhood Graph")
    url = f"{BASE_URL}/neighbors"

    # Payload matching RunNeighborsRequest DTO
    payload = {
        "project_id": TEST_PROJECT_ID,
        "dataset_id": TEST_DATASET_ID,
        "n_pcs": 30,  # Number of PCs to use (chosen from Elbow Plot)
        "n_neighbors": 3,  # Number of neighbors to compute
        "metric": "euclidean"  # Distance metric
    }

    print(f"POST {url}")
    print("Payload Config:", json.dumps(payload, indent=2))

    try:
        response = requests.post(url, json=payload)

        # We expect 201 Created
        if response.status_code == 201:
            data = response.json()
            print("\n[Success] Neighborhood Graph Built!")
            print(f"Snapshot ID: {data['snapshot_id']}")
            print(f"Saved Path: {data['snapshot_path']}")

        elif response.status_code == 404:
            print(f"\n[Failed] 404 Not Found")
            print("Server Message:", response.json().get('detail'))
            print("Tip: Did you run PCA analysis first to generate the required snapshot?")

        elif response.status_code == 500:
            print(f"\n[Failed] 500 Internal Error")
            print("Error Details:", response.text)

        else:
            print(f"\n[Failed] Status Code: {response.status_code}")
            print("Response:", response.text)

    except requests.exceptions.ConnectionError:
        print(f"\n[Error] Could not connect to backend at {url}")
        print("Is the server running on port 8888?")


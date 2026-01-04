import requests
import json
import os

API_URL = "http://127.0.0.1:8888/api/v1/projects/import/"
# !!! Modify to your real raw anndata absolute path in operating system
LOCAL_FILE_PATH = r"D:\OneDrive - Brown University\Learning\Y6\Bio Agent\datasets\test1.h5ad"


def test_simulate_user_action():
    """
    Simulates a user filling out the 'New Project' form and clicking 'Create'.
    """

    # 1. Check if file exists (Frontend validation)
    if not os.path.exists(LOCAL_FILE_PATH):
        print(f"[Error] File not found: {LOCAL_FILE_PATH}")
        return

    # 2. Prepare the JSON Payload
    # This matches the 'CreateProjectRequest' Pydantic model in the backend
    payload = {
        "project_name": "PBMC 3k Analysis Demo",
        "local_file_path": LOCAL_FILE_PATH,
        "organism": "Human",
        "tissue_type": "PBMC",
        "description": "Simulation of frontend request via Python script."
    }

    print(f"--- Sending POST Request to {API_URL} ---")
    print(f"Payload: {json.dumps(payload, indent=2)}")

    try:
        # 3. Send the Network Request
        response = requests.post(API_URL, json=payload)

        # 4. Handle Response
        if response.status_code == 201:
            print("\n[Success] Project Created Successfully!")
            print("Response Data:")
            print(json.dumps(response.json(), indent=4))
        else:
            print(f"\n[Failed] Status Code: {response.status_code}")
            print("Error Details:", response.text)

    except requests.exceptions.ConnectionError:
        print(f"\n[Error] Could not connect to backend at {API_URL}")
        print("Is the server running? (Check if uvicorn is active on port 8888)")
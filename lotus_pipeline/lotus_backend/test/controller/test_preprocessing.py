import requests
import json

# 1. Define the API Endpoint
QC_API_URL = "http://127.0.0.1:8888/api/v1/preprocessing/qc/calculate"


def test_simulate_qc_request():
    """
    Simulates a user clicking 'Run QC' in the frontend.
    Sends a POST request to trigger the QC calculation service.
    """

    TEST_PROJECT_ID = "p_1767495922_c912e"
    TEST_DATASET_ID = "ds_20260103_220522_0f5c3"

    # Prepare the JSON Payload
    # This must match the 'CalculateQCRequest' Pydantic model
    payload = {
        "project_id": TEST_PROJECT_ID,
        "dataset_id": TEST_DATASET_ID,
        "organism": "Human",
    }

    print(f"--- Sending QC Request to {QC_API_URL} ---")
    print(f"Payload: {json.dumps(payload, indent=2)}")

    try:
        # 4. Send the Network Request
        response = requests.post(QC_API_URL, json=payload)

        # 5. Handle Response
        if response.status_code == 200:
            print("\n[Success] QC Calculation Completed!")
            print("Response Data (QCResultDTO):")

            # This should print the paths to the generated JSON and PDFs
            response_data = response.json()
            print(json.dumps(response_data, indent=4))

            # Validation tip
            print("\n[Tip] You can verify these files exist in your workspace:")
            print(f" - Metrics: {response_data.get('metrics_json_path')}")
            print(f" - Violin:  {response_data.get('violin_plot_path')}")

        else:
            print(f"\n[Failed] Status Code: {response.status_code}")
            print("Error Details:", response.text)

    except requests.exceptions.ConnectionError:
        print(f"\n[Error] Could not connect to backend at {QC_API_URL}")
        print("Is the server running on port 8888?")


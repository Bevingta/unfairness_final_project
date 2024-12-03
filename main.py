import requests
import json
import os
import numpy as np
import matplotlib.pyplot as pl

# file_endpt = 'https://api.gdc.cancer.gov/files/'
# file_uuid = 'cb92f61d-041c-4424-a3e9-891b7545f351'
# response = requests.get(file_endpt + file_uuid)

# # OUTPUT METHOD 1: Write to a file.
# file = open("sample_request.json", "w")
# file.write(response.text)
# file.close()

# # OUTPUT METHOD 2: View on screen.
# print(json.dumps(response.json(), indent=2))

# Base API URL
url = 'https://api.gdc.cancer.gov/files'

# Filters & Parameters
filters = {
    "op": "or",
    "content": [
        {"op": "in",
         "content": {
             "field": "cases.project.project_id",
             "value": ["TCGA-LUAD", "TCGA-LUSC"] # Lung Adenocarcinoma (LUAD) & Lung Squamous Cell Carcinoma (LUSC)
             }
         },
        {"op": "in",
         "content": {
             "field": "cases.project.project_id",
             "value": ["TCGA-SKCM"] # Skin Cutaneous Melanoma (SKCM)
             }
         },
    ]
}

params = {
    "filters": filters,
    "fields": "id, file_name, cases.project.project_id, data_category, data_type",
    "format": "JSON",
    "size": "100" # Increase for more results
}

# Sends a request to GCD API, stored in data
response = requests.post(url, json=params)
data = response.json()

# DEBUG
# for file in data['data']['hits']:
#     print(f"ID: {file['id']}, Name: {file['file_name']}, Size: {file['file_size']}, State: {file['state']}")

# Holds file ids
file_ids = [file["id"] for file in data["data"]["hits"]]

# DEBUG
# print(file_ids)

os.makedirs("Downloads", exist_ok=True)

for id in file_ids:
    download_url = f"https://api.gdc.cancer.gov/data/{id}" # /data/ used to download data per API documentation
    response = requests.get(download_url, stream=True)

    file_name = f"{id}.tar.gz"

    # Checks if file is already downloaded
    downloaded = open("downloaded.txt", "a+", encoding="utf-8")
    downloaded.seek(0, 0) # Moves file pointer to beginning for reading
    downloaded_files = downloaded.readlines() # Stores all downloaded file names

    known_cases = []

    for file in downloaded_files:
        known_cases.append(file.strip())

    # Moves file pointer to end for appending
    downloaded.seek(0, 2)

    file_path = os.path.join("Downloads", file_name)
    with open(file_path, 'wb') as f:
        if file_name not in known_cases:
            downloaded.write(file_name + "\n")

            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)

    print(f"Downloaded: {file_name}")

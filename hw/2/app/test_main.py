from fastapi.testclient import TestClient
import pytest 
import requests
import time
from main import app, molecules_db
from app import lambda_handler
from models import Molecules
import io

client = TestClient(app)


@pytest.fixture(autouse=True)
def reset_db():
    molecules_db.clear()
    molecules_db.update({
        "m1": Molecules(identifier="m1", smile="CCO"),
        "m2": Molecules(identifier="m2", smile="CCN"),
        "m3": Molecules(identifier="m3", smile="CCCO"),
        "m4": Molecules(identifier="m4", smile="CC(=O)O"),
        "m5": Molecules(identifier="m5", smile="C1=CC=CC=C1")
    })


def test_web1_server_id():
    response = requests.get("http://localhost:8000/", headers={"Cache-Control": "no-cache"})
    assert response.status_code == 200
    assert response.json() == {"server_id": "SERVER-1"}

def test_web2_server_id():
    time.sleep(2)
    response = requests.get("http://localhost:8000/", headers={"Cache-Control": "no-cache"})
    assert response.status_code == 200
    assert response.json() == {"server_id": "SERVER-2"}


def test_retrieve_molecules():
    response = client.get("/smiles")
    assert response.status_code == 200
    assert len(response.json()) == 5


def test_retrieve_molecule_existing():
    response = client.get("/smiles/m1")
    assert response.status_code == 200
    assert response.json() == {"identifier": "m1", "smile": "CCO"}


def test_retrieve_molecule_non_existing():
    response = client.get("/smiles/non_existing")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_add_molecule():
    molecule_data = {"identifier": "m6", "smile": "CCNC"}
    response = client.post("/add", json=molecule_data)
    assert response.status_code == 200
    assert response.json() == molecule_data


def test_add_existing_molecule():
    molecule_data = {"identifier": "m1", "smile": "CCO"}
    response = client.post("/add", json=molecule_data)
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule with this identifier already exists"}


def test_update_molecule():
    updated_molecule_data = {"identifier": "m1", "smile": "CCN"}
    response = client.put("/smiles/m1", json=updated_molecule_data)
    assert response.status_code == 200
    assert response.json() == updated_molecule_data


def test_update_non_existing_molecule():
    updated_molecule_data = {"identifier": "m10", "smile": "CCN"}
    response = client.put("/smiles/m10", json=updated_molecule_data)
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_delete_molecule():
    response = client.delete("/delete/m2")
    assert response.status_code == 200
    assert response.json() == {"detail": "Molecule deleted"}


def test_delete_non_existing_molecule():
    response = client.delete("/delete/non_existing")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_substructure_search():
    response = client.get("/substructures")
    assert response.status_code == 200
    expected_result = [
        {
        "identifier": "m3",
        "substructures": [
            {
                "identifier": "m1",
                "smile": "CCO"
            }
        ]
    },
    {
        "identifier": "m4",
        "substructures": [
            {
                "identifier": "m1",
                "smile": "CCO"
            }
        ]
    }
    ]
    assert response.json() == expected_result


def test_upload_valid_csv():
    csv_content = """identifier\tsmile
    m6\tCCC
    m7\tCCOCC
    """
    files = {"file": ("test.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    response = client.post("/upload", files=files)
    assert response.status_code == 200
    assert len(response.json()["added_molecules"]) == 2
    assert "m6" in molecules_db
    assert "m7" in molecules_db


def test_upload_invalid_file_format():
    response = client.post("/upload", files={"file": ("test.txt", io.BytesIO(b"some text content"), "text/plain")})
    assert response.status_code == 400
    assert response.json() == {"detail": "Only CSV files are supported"}


def test_upload_csv_missing_columns():
    csv_content = """id\tsmile
    m6\tCCC
    """
    files = {"file": ("test.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    response = client.post("/upload", files=files)
    assert response.status_code == 400
    assert response.json() == {"detail": "File must contain 'identifier' and 'smile' columns"}


def test_upload_csv_duplicate_entries():
    csv_content = """identifier\tsmile
    m1\tCCC
    """
    files = {"file": ("test.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    response = client.post("/upload", files=files)
    assert response.status_code == 200
    assert len(response.json()["added_molecules"]) == 0
    assert molecules_db["m1"].smile == "CCO"  # Ensure original entry is unchanged


def test_upload_csv_with_invalid_smiles():
    csv_content = """identifier\tsmile
    m8\tINVALID
    """
    files = {"file": ("test.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    response = client.post("/upload", files=files)
    assert response.status_code == 200
    assert len(response.json()["added_molecules"]) == 1
    assert "m8" in molecules_db
    assert molecules_db["m8"].smile == "INVALID"  
    

def test_lambda_handler():
    event = {"name": "Test"}
    result = lambda_handler(event, None)
    assert result['statusCode'] == 200
    assert result['body'] == 'Hello, Test!'

def test_lambda_handler_default():
    event = {}
    result = lambda_handler(event, None)
    assert result['statusCode'] == 200
    assert result['body'] == 'Hello, World!'
from fastapi import FastAPI, HTTPException, UploadFile, File
from models import Molecules
from rdkit import Chem
import io
import os
import csv

app = FastAPI()

molecules_db = {
    "m1": Molecules(identifier="m1", smile="CCO"),
    "m2": Molecules(identifier="m2", smile="CCN"),
    "m3": Molecules(identifier="m3", smile="CCCO"),
    "m4": Molecules(identifier="m4", smile="CC(=O)O"),
    "m5": Molecules(identifier="m5", smile="C1=CC=CC=C1")
}

@app.get("/")
def get_server():
    return {"server_id": os.getenv("SERVER_ID", "1")}


@app.get("/smiles")
def retrive_molecules():
    """
    Endpoint to retrieve all smiles from the database.
    Returns a dictionary with all smiles.
    """
    return molecules_db


@app.get("/smiles/{identifier}")
def retrieve_molecule(identifier: str):
    """
    Endpoint to retrieve a specific smile by its identifier.
    Args:
        identifier: The unique identifier of the smile.
    Returns:
        The smile with the given identifier if it exists.
    Raises:
        HTTPException: If the smile with the given identifier does not exist.
    """
    molecule = molecules_db.get(identifier)
    if molecule:
        return molecule
    raise HTTPException(status_code=404, detail="Molecule not found")


@app.post("/add")
def add_molecule(molecule: Molecules):
    """
    Endpoint to add a new smile to the database.
    Args:
        molecule: The smile object to be added.
    Returns:
        The added smile.
    Raises:
        HTTPException: If a smile with the same identifier already exists.
    """
    if molecule.identifier in molecules_db:
        raise HTTPException(status_code=400, detail="Molecule with this identifier already exists")
    molecules_db[molecule.identifier] = molecule
    return molecule


@app.put("/smiles/{identifier}")
def update_molecule(identifier: str, updated_molecule: Molecules):
    """
    Endpoint to update an existing smile's information.
    Args:
        identifier: The unique identifier of the smile to be updated.
        updated_molecule: The updated smile object.
    Returns:
        The updated smile.
    Raises:
        HTTPException: If the smile to be updated does not exist or if the new identifier already exists.
    """
    if identifier in molecules_db:
        # Remove the old entry if the identifier is different
        if updated_molecule.identifier != identifier:
            if updated_molecule.identifier in molecules_db:
                raise HTTPException(status_code=400, detail="New identifier already exists")
            del molecules_db[identifier]
            molecules_db[updated_molecule.identifier] = updated_molecule
        else:
            molecules_db[identifier] = updated_molecule
        return updated_molecule
    raise HTTPException(status_code=404, detail="Molecule not found")



@app.delete("/delete/{identifier}")
def delete_molecule(identifier: str):
    """
    Endpoint to delete a specific smile by its identifier.
    Args:
        identifier: The unique identifier of the smile to be deleted.
    Returns:
        A confirmation message upon successful deletion.
    Raises:
        HTTPException: If the smile with the given identifier does not exist.
    """
    if identifier in molecules_db:
        del molecules_db[identifier]
        return {"detail": "Molecule deleted"}
    raise HTTPException(status_code=404, detail="Molecule not found")


@app.get("/substructures")
def substructure_search():
    """
    Search for substructures within the molecules database.
    Returns a list of molecules that contain substructures from other molecules.
    """
    # Create a dictionary to store molecule identifiers and their molecules
    molecules_dict = {mol.identifier: mol for mol in molecules_db.values()}

    results = []
    substructure_matches = {identifier: [] for identifier in molecules_dict}

    # Create molecule objects for all molecules once
    molecule_objects = {identifier: Chem.MolFromSmiles(molecule.smile) 
                         for identifier, molecule in molecules_dict.items()}

    # Check substructure matches
    for identifier, mol in molecule_objects.items():
        if mol is None:
            continue

        for sub_id, sub_mol in molecule_objects.items():
            if identifier == sub_id or sub_mol is None:
                continue

            if mol.HasSubstructMatch(sub_mol):
                substructure_matches[identifier].append({
                    "identifier": sub_id,
                    "smile": molecules_dict[sub_id].smile
                })

    # Prepare the final result list
    for identifier, matches in substructure_matches.items():
        if matches:
            results.append({
                "identifier": identifier,
                "substructures": matches
            })

    return results


@app.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    """
    Endpoint to upload a CSV file containing smile data.
    The CSV file must have 'identifier' and 'smile' columns and use a tab delimiter.
    Args:
        file: The uploaded CSV file.
    Returns:
        A list of added smile.
    Raises:
        HTTPException: If the file format is incorrect or if there are issues with the CSV content.
    """
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="Only CSV files are supported")
    
    content = await file.read()
    data = io.StringIO(content.decode("utf-8"))

    try:
        csv_reader = csv.DictReader(data, delimiter='\t')
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to read CSV: {str(e)}")

    added_smiles = []

    for row in csv_reader:
        identifier, smile = row.get('identifier'), row.get('smile')
        if not identifier or not smile:
            raise HTTPException(status_code=400, detail="File must contain 'identifier' and 'smile' columns")

        if identifier not in molecules_db:
            molecule = Molecules(identifier=identifier, smile=smile)
            molecules_db[identifier] = molecule
            added_smiles.append(molecule)

    return {"added_molecules": added_smiles}
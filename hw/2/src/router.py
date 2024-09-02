from fastapi import APIRouter, HTTPException, UploadFile, File, Depends
from config import SessionLocal
from sqlalchemy.orm import Session
from schemas import MoleculeSchema, RequestMolecule, Response
import crud
from logger import logger
from typing import List
from substructure_search import substructure

router = APIRouter()


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@router.get("/smiles", response_model=Response[List[MoleculeSchema]])
async def get(db: Session = Depends(get_db)):
    """
    Endpoint to retrieves all smiles from the database.
    """
    molecules = crud.get_molecules(db)
    
    molecules_list = list(molecules)  # Convert generator to a list for logging

    logger.info(f"Retrieved {len(molecules_list)} molecules.")

    return Response(code="200", status="Success", message="Molecules retrieved", result=molecules_list)


@router.get("/smiles/{offset}", response_model=Response[List[MoleculeSchema]])
async def get(db: Session = Depends(get_db), limit: int = 5, offset: int = 0):
    """
    Endpoint to retrieves 5 smiles from the database.
    Returns smiles based on page each page has 5 smiles, first page is 0.
    """
    offset = limit * offset
    
    molecules = crud.get_molecules(db, limit, offset)
    
    molecules_list = list(molecules)  # Convert generator to a list for logging

    logger.info(f"Retrieved {len(molecules_list)} molecules with limit={limit}")

    return Response(code="200", status="Success", message="Molecules retrieved", result=molecules_list)


@router.get("/smile/{identifier}", response_model=Response[MoleculeSchema])
async def get(identifier: str, db: Session = Depends(get_db)):
    """
    Endpoint to retrieve a specific smile by its identifier.
    Args:
        identifier: The unique identifier of the smile.
    Returns:
        The smile with the given identifier if it exists.
    Raises:
        HTTPException: If the smile with the given identifier does not exist.
    """
    molecule = crud.get_molecule_by_id(db, identifier)
    if molecule:
        return Response(code="200", status="Success", message="Molecule retrieved", result=molecule)
    raise HTTPException(status_code=404, detail="Molecule not found")


@router.put("/smile/{identifier}", response_model=Response[MoleculeSchema])
async def put(identifier: str, molecule: RequestMolecule, db: Session = Depends(get_db)):
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
    try:
        updated_molecule = crud.update_molecule(db, identifier, molecule.parameters)
        
        if updated_molecule:
            return Response(code="200", status="Success", message="Molecule updated", result=updated_molecule)
        
        raise HTTPException(status_code=404, detail="Molecule not found")
    
    except ValueError:
        raise HTTPException(status_code=400, detail="Molecule with identifier already exists")


@router.delete("/smile/{identifier}", response_model=Response[None])
async def delete_molecule(identifier: str, db: Session = Depends(get_db)):
    """
    Endpoint to delete a specific molecule by its identifier.

    Args:
        identifier (str): The unique identifier of the molecule to be deleted.

    Returns:
        Response[None]: A confirmation message upon successful deletion.

    Raises:
        HTTPException: If the molecule with the given identifier does not exist.
    """
    # Attempt to delete the molecule
    result = crud.delete_molecule(db, identifier)
    
    # Check if the deletion was successful
    if result:
        return Response(code="200", status="Success", message="Molecule successfully deleted", result=None)
    
    # If the molecule was not found, raise a 404 error
    raise HTTPException(status_code=404, detail="Molecule not found")


@router.post("/add", response_model=Response[MoleculeSchema])
async def post(molecule: RequestMolecule, db: Session = Depends(get_db)):
    """
    Endpoint to add a new smile to the database.
    Args:
        molecule: The smile object to be added.
    Returns:
        The added smile.
    Raises:
        HTTPException: If a smile with the same identifier already exists.
    """
    existing_molecule = crud.get_molecule_by_id(db, molecule.parameters.identifier)
    if existing_molecule:
        raise HTTPException(status_code=400, detail="Molecule with this identifier already exists.")
    
    new_molecule = crud.add_molecule(db, molecule.parameters)
    return Response(code="201", status="Success", message="Molecule added", result=new_molecule)


@router.get("/substructures")
def substructure_search(db: Session = Depends(get_db)):
    """
    Search for substructures within the molecules database.
    Returns a list of molecules that contain substructures from other molecules.
    """
    return substructure(db)


@router.post("/add/csv", response_model=Response[None])
async def post_csv(csv_file: UploadFile = File(...), db: Session = Depends(get_db)):
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
    try:
        result = crud.add_molecule_from_csv(db, csv_file)
        return Response(code="201", status="Success", message="Molecules added from CSV", result=result)
    except HTTPException as e:
        raise HTTPException(status_code=400, detail=f"Unknown: {e}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Not Valid: {e}")

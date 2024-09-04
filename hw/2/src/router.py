from fastapi import APIRouter, HTTPException, UploadFile, File, Depends
from config import SessionLocal
from sqlalchemy.orm import Session
from schemas import MoleculeSchema, RequestMolecule, Response
import crud
import redis
import json
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


# Initialized Redis client
redis_client = redis.Redis(host='redis', port=6379, db=0)

# Set the cache expiration time to 1 hour
CACHE_EXPIRATION_TIME = 3600  # 1 hour

def get_cached_result(key: str):
    '''Get the cached result from Redis'''
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: dict, expiration: int = 60):
    '''Set the cache in Redis'''
    redis_client.setex(key, expiration, json.dumps(value))


def fetch_and_cache_data(cache_key: str, fetch_function, *args):
    """Fetch data from DB, cache it, and return."""
    # Check if cached data exists
    cached_data = get_cached_result(cache_key)
    if cached_data:
        logger.info(f"Data returned from cache for key: {cache_key}.")
        return cached_data, "cache"

    # Fetch data from DB
    data = fetch_function(*args)
    data_list = list(data)

    # Convert data to dictionary format
    data_dict_list = [MoleculeSchema.from_orm(item).dict() for item in data_list]

    # Cache the data
    set_cache(cache_key, data_dict_list)

    logger.info(f"Data fetched from database for key: {cache_key}.")
    return data_dict_list, "database"


@router.get("/smiles", response_model=Response[List[MoleculeSchema]])
async def get(db: Session = Depends(get_db)):
    """
    Endpoint to retrieve all smiles from the database.
    Returns a dictionary with all smiles.
    """
    cache_key = "all_smiles"
    
    # Fetch and cache data
    smiles_list, source = fetch_and_cache_data(cache_key, crud.get_molecules, db)

    return Response(code="200", status="Success", message=f"Molecules retrieved from {source}", result=smiles_list)


@router.get("/smiles/{offset}", response_model=Response[List[MoleculeSchema]])
async def get(db: Session = Depends(get_db), limit: int = 5, offset: int = 0):
    """
    Endpoint to retrieves 5 smiles from the database.
    Returns smiles based on page each page has 5 smiles, first page is 0.
    """
    offset = limit * offset
    cache_key = f"smiles_{offset}_{limit}"
    
    # Fetch and cache data
    smiles_list, source = fetch_and_cache_data(cache_key, crud.get_molecules, db, limit, offset)

    return Response(code="200", status="Success", message=f"Molecules retrieved from {source}", result=smiles_list)



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

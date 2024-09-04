from fastapi import HTTPException
from sqlalchemy.orm import Session
from models import Molecule
from schemas import MoleculeSchema
from logger import logger
import io
import csv


def get_molecules(db: Session, limit: int = None, offset: int = 0):
    logger.info("Fetching all molecules from the database")
    
    if limit:
        query = db.query(Molecule).offset(offset).limit(limit)
    else:
        query = db.query(Molecule)
        
    for molecule in query:
        logger.debug(f"Found molecule with identifier {molecule.identifier}")
        yield molecule
    return query

def get_molecule_by_id(db: Session, identifier: str):
    logger.info(f"Fetching molecule with identifier {identifier} from the database")
    query = db.query(Molecule).filter(Molecule.identifier == identifier).first()
    if query:
        logger.debug(f"Found molecule with identifier {identifier}")
        return query
    logger.warning(f"Molecule with identifier {identifier} not found")
    return None

def update_molecule(db: Session, identifier: str, updated_molecule: MoleculeSchema):
    logger.info(f"Updating molecule with identifier {identifier}")
    molecule = db.query(Molecule).filter(Molecule.identifier == identifier).first()
    
    if not molecule:
        logger.warning(f"Molecule with identifier {identifier} not found")
        return None
    
    if updated_molecule.identifier != identifier:
        existing_molecule = db.query(Molecule).filter(Molecule.identifier == updated_molecule.identifier).first()
        if existing_molecule:
            logger.error(f"Molecule with identifier {updated_molecule.identifier} already exists")
            raise ValueError("A molecule with the new identifier already exists.")
    
    molecule.identifier = updated_molecule.identifier
    molecule.smile = updated_molecule.smile
    db.commit()
    logger.debug(f"Molecule with identifier {identifier} updated successfully")
    return molecule

def delete_molecule(db: Session, identifier: str):
    logger.info(f"Deleting molecule with identifier {identifier}")
    molecule = db.query(Molecule).filter(Molecule.identifier == identifier).first()
    if molecule:
        db.delete(molecule)
        db.commit()
        logger
        return {"detail": "Molecule deleted"}
    logger.warning(f"Molecule with identifier {identifier} not found")
    return None

def add_molecule(db: Session, molecule: MoleculeSchema):
    logger.info(f"Adding molecule with identifier {molecule.identifier}")
    new_molecule = Molecule(identifier=molecule.identifier, smile=molecule.smile)
    if new_molecule:
        db.add(new_molecule)
        db.commit()
        logger.debug(f"Molecule with identifier {molecule.identifier} added successfully")
        return new_molecule
    logger.error(f"Failed to add molecule with identifier {molecule.identifier}")
    return None

def add_molecule_from_csv(db: Session, csv_file):
    logger.info("Adding molecules from CSV file")
    # Ensure the file has a .csv extension
    if not csv_file.filename.endswith('.csv'):
        logger.error("Only CSV files are supported")
        raise HTTPException(status_code=400, detail="Only CSV files are supported")

    # Read the file content
    content = csv_file.file.read()
    data = io.StringIO(content.decode("utf-8"))

    try:
        # Parse the CSV data
        csv_reader = csv.DictReader(data, delimiter='\t')
    except Exception as e:
        logger.error(f"Failed to read CSV: {str(e)}")
        raise HTTPException(status_code=400, detail=f"Failed to read CSV: {str(e)}")

    added_smiles = []

    for row in csv_reader:
        identifier, smile = row.get('identifier'), row.get('smile')

        # Check for missing data
        if not identifier or not smile:
            logger.error("File must contain 'identifier' and 'smile' columns")  
            raise HTTPException(status_code=400, detail="File must contain 'identifier' and 'smile' columns")

        # Check if the molecule already exists
        existing_molecule = db.query(Molecule).filter(Molecule.identifier == identifier).first()
        if not existing_molecule:
            logger.debug(f"Adding molecule with identifier {identifier}")
            new_molecule = Molecule(identifier=identifier, smile=smile)
            added_smiles.append(new_molecule)
            db.add(new_molecule)

    # Commit all the changes
    db.commit()

    # Return None if no molecules were added
    if not added_smiles:
        logger.warning("No molecules were added")
        return None

    # Return the added molecules
    return {"added_molecules": added_smiles}
    
    # try:
    #     csv_content = io.StringIO(csv_file.file.read().decode("utf-8"))
    #     csv_reader = csv.DictReader(csv_content, delimiter='\t')
    # except Exception as e:
    #     raise HTTPException(status_code=400, detail=f"Failed to read CSV: {str(e)}")
    
    # added_molecules = []

    # # Process each row in the CSV file
    # for row in csv_reader:
    #     identifier = row.get('identifier')
    #     smile = row.get('smile')

    #     # Validate that both identifier and smile columns exist
    #     if not identifier or not smile:
    #         raise HTTPException(status_code=400, detail="CSV must contain 'identifier' and 'smile' columns")

    #     # Check if the molecule with the same identifier already exists
    #     existing_molecule = db.query(Molecule).filter(Molecule.identifier == identifier).first()
    #     if not existing_molecule:
    #         new_molecule = Molecule(identifier=identifier, smile=smile)
    #         db.add(new_molecule)
    #         added_molecules.append(new_molecule)
    
    # # Commit all the changes to the database at once
    # db.commit()

    # return {"added_molecules": added_molecules}

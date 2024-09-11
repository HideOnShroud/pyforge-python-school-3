from sqlalchemy.orm import Session
from rdkit import Chem
import crud

def substructure(db: Session):
    # Fetch all molecules from the database
    molecules = crud.get_molecules(db)
    
    if not molecules:
        return []

    # Convert the list of molecule objects to a dictionary with identifiers as keys
    molecules_dict = {mol.identifier: mol for mol in molecules}

    # Create RDKit molecule objects for all molecules
    def create_molecule_objects():
        for identifier, mol in molecules_dict.items():
            try:
                yield identifier, Chem.MolFromSmiles(mol.smile)
            except Exception as e:
                # Log or handle the invalid SMILES string here
                print(f"Error processing SMILES for {identifier}: {e}")
                yield identifier, None

    molecule_objects = dict(create_molecule_objects())

    substructure_matches = {identifier: [] for identifier in molecule_objects if molecule_objects[identifier] is not None}

    def find_substructure_matches():
        for identifier, mol in molecule_objects.items():
            if mol is None:
                continue

            for sub_id, sub_mol in molecule_objects.items():
                if identifier == sub_id or sub_mol is None:
                    continue

                if mol.HasSubstructMatch(sub_mol):
                    yield identifier, {
                        "identifier": sub_id,
                        "smile": molecules_dict[sub_id].smile
                    }

    for identifier, match in find_substructure_matches():
        substructure_matches[identifier].append(match)

    results = [
        {
            "identifier": identifier,
            "substructures": matches
        }
        for identifier, matches in substructure_matches.items() if matches
    ]

    return results

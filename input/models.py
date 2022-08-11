from fastapi import HTTPException
from pydantic import BaseModel, validator
from rdkit import Chem
from rdkit.Chem import AllChem


class InputValidator(BaseModel):
    reaction_smarts: str
    reactants: str

    @validator("reaction_smarts")
    def check_if_smarts_valid(cls, value):
        AllChem.ReactionFromSmarts(value)
        return value
        
    @validator("reactants")
    def check_if_smiles_valid(cls, value):
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise HTTPException(status_code=422, detail="SMILES is not correct.")
        return value

from fastapi import HTTPException
from pydantic import BaseModel, validator
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List, Union


class InputValidator(BaseModel):
    reaction_smarts: str
    reactants: Union[str, List[str]]

    @validator("reaction_smarts")
    def check_if_smarts_valid(cls, value):
        AllChem.ReactionFromSmarts(value)
        return value
        
    @validator("reactants")
    def check_if_smiles_valid(cls, value):
        if isinstance(value, str):
            value = (value,)

        for code in value:
            molecule = Chem.MolFromSmiles(code)
            if molecule is None:
                raise HTTPException(status_code=422, detail=f"SMILES code {code} is not correct.")

        return value

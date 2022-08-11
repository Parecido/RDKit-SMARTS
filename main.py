from fastapi import FastAPI, HTTPException
from input.models import InputValidator
from rdkit import Chem
from rdkit.Chem import AllChem


def transform(smarts, smiles):
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        ps = rxn.RunReactants((Chem.MolFromSmiles(smiles),))
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
        
    return ps


def analyse_reactions(ps):
    smiles_valid = []
    smiles_invalid = []
    
    for reaction in ps:
        products = []
        
        for reactant in reaction:
            try:
                Chem.SanitizeMol(reactant)
                products.append(reactant)
            except Exception:
                break
                
        if len(products) == len(reaction):
            smiles_valid.append(Chem.MolToSmiles(reactant) for reactant in reaction)
        else:
            smiles_invalid.append(Chem.MolToSmiles(reactant) for reactant in reaction)

    unique_smiles_valid = [code for code in set(tuple(code) for code in smiles_valid)]
    unique_smiles_invalid = [code for code in set(tuple(code) for code in smiles_invalid)]

    return unique_smiles_valid, unique_smiles_invalid


app = FastAPI()


@app.get("/smarts/")
async def run_smarts(request: InputValidator):
    reaction_smarts = request.reaction_smarts
    reactants = request.reactants
    
    ps = transform(reaction_smarts, reactants)
    valid, invalid = analyse_reactions(ps)
    
    response = {
        "products_valid": valid,
        "products_invalid": invalid,
    }

    return response

from fastapi import FastAPI, HTTPException
from input.models import InputValidator
from rdkit import Chem
from rdkit.Chem import AllChem


def transform(smarts, smiles):
    try:
        if isinstance(smiles, str):
            smiles = (smiles,)

        mols = ([Chem.MolFromSmiles(code) for code in smiles])
        rxn = AllChem.ReactionFromSmarts(smarts)
        ps = rxn.RunReactants(mols)
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
        
    return ps


def analyse_reactions(ps):
    smiles_valid = []
    smiles_invalid = []
    
    for reaction in ps:
        products = []
        
        for product in reaction:
            try:
                Chem.SanitizeMol(product)
                products.append(product)
            except Exception:
                break
                
        if len(products) == len(reaction):
            smiles_valid.append(Chem.MolToSmiles(product) for product in reaction)
        else:
            smiles_invalid.append(Chem.MolToSmiles(product) for product in reaction)

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

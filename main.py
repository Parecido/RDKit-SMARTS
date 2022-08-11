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
    smiles = []
    
    for reaction in ps:
        products = []
        
        for reactant in reaction:
            try:
                Chem.SanitizeMol(reactant)
                products.append(Chem.MolToSmiles(reactant))
            except ValueError:
                break
                
        if len(products) == len(reaction):
            smiles.append(products)
        
    unique_smiles = [code for code in set(tuple(code) for code in smiles)]
    return unique_smiles


app = FastAPI()


@app.get("/smarts/")
async def run_smarts(request: InputValidator):
    reaction_smarts = request.reaction_smarts
    reactants = request.reactants
    
    ps = transform(reaction_smarts, reactants)
    smiles = analyse_reactions(ps)
    
    response = { 
        "products_smiles": smiles,
    }

    return response

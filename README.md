# RDKit-SMARTS

## Getting Started

### Dependencies

* Python 3.8+ (Anaconda)
* OS: Linux, MacOS
* Libraries: pydantic, RDKit, FastAPI

### Installing

* Download latest RDKit-SMARTS version
* Prepare conda environment
```
conda create -n rdkit-smarts --file=requirements.txt
```

### Executing program

#### Executing RDKit-SMARTS server

```
conda activate rdkit-smarts
uvicorn main:app --host 127.0.0.1 --port 8000
```

#### Executing RDKit-SMARTS client

```
import requests

data = {
    'reaction_smarts': '[c:8]-[c:6]>>[c:8][I:55].[B:99][c:6]',
    'reactants': 'CC1=CC=C(C=C1)C1=CC(=CC=C1C)C1=CC(C)=CC(C)=C1'
}
    
r = requests.get("http://127.0.0.1:8000/smarts/", json=data)
print(r.json())
```

### Data structures

#### Input structure

The input file should contain the "reaction_smarts" and "reactants" keys, which correspond to the SMARTS reaction and SMILES molecular code, respectively. Both keys are of type string.
```
{
    'reaction_smarts': '[c:8]-[c:6]>>[c:8][I:55].[B:99][c:6]',
    'reactants': 'CC1=CC=C(C=C1)C1=CC(=CC=C1C)C1=CC(C)=CC(C)=C1'
}
```

#### Output structure

The RDKit-SMARTS program sends json with "products_valid" and "products_invalid" keys with SMILES code. For more details see "Wrong structures" section.
```
{
    'products_valid': [['Cc1cc(C)cc(-c2ccc(C)c(I)c2)c1', 'Bc1ccc(C)cc1'], ['Cc1ccc(I)cc1', 'Bc1cc(-c2cc(C)cc(C)c2)ccc1C'], ['Cc1cc(C)cc(I)c1', 'Bc1ccc(C)c(-c2ccc(C)cc2)c1'], ['Cc1ccc(-c2cc(I)ccc2C)cc1', 'Bc1cc(C)cc(C)c1']],
    'products_invalid': []
}
```

### Output limits

#### Molecular symmetry

SMARTS can generate several identical solutions for the same molecule. This is due to the molecular symmetry, which modifies the molecule in several identical ways. Therefore, RDKit-SMARTS only generates a list of unique products.

#### Wrong structures

SMARTS can sometimes lead to a chemically questionable product structure. To avoid this problem, sanitization has been applied. If RDKit-SMARTS obtains a structure that is not chemically represented and cannot be corrected, then it will not be sent to the output.

## Authors

Michal Michalski ([Git](https://github.com/Parecido), [LinkedIn](https://www.linkedin.com/in/michal-michalski95), [Orcid](https://orcid.org/0000-0001-6969-2074), [ResearchGate](https://www.researchgate.net/profile/Michal-Michalski-5))

## TODO

* Support more than one reagent
* Stereochemistry (RDKit does not support non-tetrahedral chiral classes)
* Test cis-trans isomerism (The new version of RDKit should supports it)


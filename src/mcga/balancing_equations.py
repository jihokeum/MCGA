##to balance equations
"""on bash pip install chempy"""

from chempy import balance_stoichiometry
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def smiles_to_formula(smiles):
    """Convert SMILES to a molecular formula (Hill system) using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return rdMolDescriptors.CalcMolFormula(mol)

def get_balanced_equation(reactant_smiles, product_smiles):
    """
    Balances a chemical reaction from SMILES and returns coefficients and formulas.

    Returns:
        dict with keys:
        - 'reactants': {formula: coefficient}
        - 'products' : {formula: coefficient}
        - 'formatted': human-readable equation
        - 'formula_to_smiles': {formula: SMILES}
    """
    reactant_formulas = {}
    product_formulas = {}
    smiles_to_formula_map = {}

    for s in reactant_smiles:
        f = smiles_to_formula(s)
        if f:
            reactant_formulas[f] = reactant_formulas.get(f, 0) + 1
            smiles_to_formula_map[f] = s

    for s in product_smiles:
        f = smiles_to_formula(s)
        if f:
            product_formulas[f] = product_formulas.get(f, 0) + 1
            smiles_to_formula_map[f] = s

    try:
        reac_bal, prod_bal = balance_stoichiometry(reactant_formulas, product_formulas)

        formatted = ' + '.join(f"{v} {k}" for k, v in reac_bal.items())
        formatted += " → "
        formatted += ' + '.join(f"{v} {k}" for k, v in prod_bal.items())

        return {
            "reactants": dict(reac_bal),
            "products": dict(prod_bal),
            "formatted": formatted,
            "formula_to_smiles": smiles_to_formula_map
        }
    except Exception as e:
        return {
            "reactants": {},
            "products": {},
            "formatted": f"Could not balance equation: {e}",
            "formula_to_smiles": {}
        }
    
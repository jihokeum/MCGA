##to balance equations
"""on bash pip install chempy"""

from chempy import balance_stoichiometry
from chempy.util.parsing import formula_to_composition

def smiles_to_formula(smiles):
    """Convert SMILES to Hill formula using RDKit."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return rdMolDescriptors.CalcMolFormula(mol)

def get_balanced_equation(reactant_smiles, product_smiles):
    """Balance a chemical equation given reactant and product SMILES strings."""
    reactant_formulas = {}
    product_formulas = {}

    for s in reactant_smiles:
        formula = smiles_to_formula(s)
        if not formula:
            continue
        reactant_formulas[formula] = reactant_formulas.get(formula, 0) + 1

    for s in product_smiles:
        formula = smiles_to_formula(s)
        if not formula:
            continue
        product_formulas[formula] = product_formulas.get(formula, 0) + 1

    try:
        reac_bal, prod_bal = balance_stoichiometry(reactant_formulas, product_formulas)
        balanced_eq = ' + '.join(f"{v} {k}" for k, v in reac_bal.items())
        balanced_eq += " → "
        balanced_eq += ' + '.join(f"{v} {k}" for k, v in prod_bal.items())
        return balanced_eq
    except Exception as e:
        return f"Could not balance equation: {e}"



import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from streamlit_ketcher import st_ketcher
from io import BytesIO

##drawing or input molecules (reactants/prodcuts/agents as drawing or smiles strings in streamlit)

st.title("Green Chemistry Reaction Input (Draw or Type)")

# Unified input: draw or type SMILES
def get_smiles_input(label, key):
    st.subheader(f"{label}")
    drawn = st_ketcher(f"Draw {label}", key=f"{key}_draw")
    typed = st.text_input(f"Or enter {label} SMILES manually:", key=f"{key}_typed")

    if drawn:
        return drawn.strip()
    elif typed:
        return typed.strip()
    else:
        return ""

# Convert Mol to image
def draw_mol(smiles, label):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        buf = BytesIO()
        img.save(buf, format="PNG")
        st.image(buf.getvalue(), caption=label)
    else:
        st.error(f"{label} SMILES is invalid or empty.")

# Get inputs for each component
reactant_smiles = get_smiles_input("Reactant", "reactant")
product_smiles = get_smiles_input("Product", "product")
agent_smiles = get_smiles_input("Agent (optional)", "agent")

# Submit and show results
if st.button("Submit", key="submit_btn"):
    reaction_data = {
        "reactants": [reactant_smiles] if reactant_smiles else [],
        "products": [product_smiles] if product_smiles else [],
        "agents": [agent_smiles] if agent_smiles else []
    }

    st.subheader("Parsed Reaction")
    st.json(reaction_data)

    # Show molecules if available
    if reactant_smiles:
        st.subheader("Reactant Structure")
        draw_mol(reactant_smiles, "Reactant")

    if product_smiles:
        st.subheader("Product Structure")
        draw_mol(product_smiles, "Product")

    if agent_smiles:
        st.subheader("Agent Structure")
        draw_mol(agent_smiles, "Agent")

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from streamlit_ketcher import st_ketcher
from io import BytesIO
import pubchempy as pcp
import requests
import time

##drawing or input molecules (reactants/prodcuts/agents as drawing or smiles strings in streamlit)

st.title("Green Chemistry Reaction Input (Draw or Type)")

# Function to convert chemical name to SMILES
def name_to_smiles(chemical_name):
    if not chemical_name:
        return None
    
    try:
        # Try with PubChemPy
        compounds = pcp.get_compounds(chemical_name, 'name')
        if compounds:
            return compounds[0].canonical_smiles
    except Exception as e:
        st.warning(f"PubChem API error: {e}")
        
    # Backup method using PubChem API directly
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/JSON"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except Exception as e:
        st.warning(f"Error in direct PubChem API: {e}")
    
    return None

# Unified input: draw, type SMILES, or enter chemical name
def get_smiles_input(label, key):
    st.subheader(f"{label}")
    
    # Create tabs for different input methods
    tab1, tab2, tab3 = st.tabs(["Draw Structure", "Enter SMILES", "Chemical Name"])
    
    with tab1:
        drawn = st_ketcher(f"Draw {label}", key=f"{key}_draw")
    
    with tab2:
        typed = st.text_input(f"Enter {label} SMILES:", key=f"{key}_typed")
    
    with tab3:
        chem_name = st.text_input(f"Enter {label} chemical name:", key=f"{key}_name")
        name_smiles = None
        if chem_name:
            # Store the chemical name for later use
            st.session_state[f"{key}_chem_name"] = chem_name
    
    # Processing priority: drawn > typed > name-derived
    if drawn:
        return drawn.strip(), "drawing"
    elif typed:
        return typed.strip(), "smiles"
    else:
        return "", "name"  # Return empty string and indicate it's from a name input

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
reactant_result = get_smiles_input("Reactant", "reactant")
reactant_smiles, reactant_source = reactant_result

product_result = get_smiles_input("Product", "product")
product_smiles, product_source = product_result

agent_result = get_smiles_input("Agent (optional)", "agent")
agent_smiles, agent_source = agent_result

# Submit and show results
if st.button("Submit", key="submit_btn"):
    with st.spinner("Processing chemical names and structures..."):
        # Process chemical names if needed
        if reactant_source == "name" and "reactant_chem_name" in st.session_state:
            reactant_name = st.session_state["reactant_chem_name"]
            reactant_smiles = name_to_smiles(reactant_name)
            if reactant_smiles:
                st.success(f"Successfully converted '{reactant_name}' to SMILES")
            else:
                st.error(f"Could not find SMILES for '{reactant_name}'")
        
        if product_source == "name" and "product_chem_name" in st.session_state:
            product_name = st.session_state["product_chem_name"]
            product_smiles = name_to_smiles(product_name)
            if product_smiles:
                st.success(f"Successfully converted '{product_name}' to SMILES")
            else:
                st.error(f"Could not find SMILES for '{product_name}'")
        
        if agent_source == "name" and "agent_chem_name" in st.session_state:
            agent_name = st.session_state["agent_chem_name"]
            agent_smiles = name_to_smiles(agent_name)
            if agent_smiles:
                st.success(f"Successfully converted '{agent_name}' to SMILES")
            else:
                st.error(f"Could not find SMILES for '{agent_name}'")
        
        # Short delay to allow the success/error messages to be displayed
        time.sleep(0.5)

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
        if reactant_source == "name" and "reactant_chem_name" in st.session_state:
            st.caption(f"Generated from chemical name: {st.session_state['reactant_chem_name']}")

    if product_smiles:
        st.subheader("Product Structure")
        draw_mol(product_smiles, "Product")
        if product_source == "name" and "product_chem_name" in st.session_state:
            st.caption(f"Generated from chemical name: {st.session_state['product_chem_name']}")

    if agent_smiles:
        st.subheader("Agent Structure")
        draw_mol(agent_smiles, "Agent")
        if agent_source == "name" and "agent_chem_name" in st.session_state:
            st.caption(f"Generated from chemical name: {st.session_state['agent_chem_name']}")
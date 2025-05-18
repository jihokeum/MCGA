import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, rdMolDescriptors
from streamlit_ketcher import st_ketcher
from io import BytesIO
import pubchempy as pcp
import requests
import time
import pandas as pd
import numpy as np

##drawing or input molecules (reactants/prodcuts/agents as drawing or smiles strings in streamlit)

st.title("Green Chemistry Reaction Input & Prediction")

# Initialize session state for dynamic components
if 'reactants_count' not in st.session_state:
    st.session_state.reactants_count = 1
if 'products_count' not in st.session_state:
    st.session_state.products_count = 1
if 'agents_count' not in st.session_state:
    st.session_state.agents_count = 1

# Function to handle adding new components
def add_component(component_type):
    if component_type == 'reactant':
        st.session_state.reactants_count += 1
    elif component_type == 'product':
        st.session_state.products_count += 1
    elif component_type == 'agent':
        st.session_state.agents_count += 1

# Function to handle removing components
def remove_component(component_type):
    if component_type == 'reactant' and st.session_state.reactants_count > 1:
        st.session_state.reactants_count -= 1
    elif component_type == 'product' and st.session_state.products_count > 1:
        st.session_state.products_count -= 1
    elif component_type == 'agent' and st.session_state.agents_count > 1:
        st.session_state.agents_count -= 1

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
def get_smiles_input(label, key, index=None):
    # Generate a unique key for each component
    unique_key = f"{key}_{index}" if index is not None else key
    
    st.subheader(f"{label} {index if index is not None else ''}")
    
    # Create tabs for different input methods
    tab1, tab2, tab3 = st.tabs(["Draw Structure", "Enter SMILES", "Chemical Name"])
    
    with tab1:
        # Utiliser une chaîne vide comme premier argument pour avoir un éditeur vide
        drawn = st_ketcher(value="", key=f"{unique_key}_draw")
    
    with tab2:
        typed = st.text_input(f"Enter {label} SMILES:", key=f"{unique_key}_typed")
    
    with tab3:
        chem_name = st.text_input(f"Enter {label} chemical name:", key=f"{unique_key}_name")
        name_smiles = None
        if chem_name:
            # Store the chemical name for later use
            st.session_state[f"{unique_key}_chem_name"] = chem_name
    
    # Processing priority: drawn > typed > name-derived
    if drawn:
        return drawn.strip(), "drawing", unique_key
    elif typed:
        return typed.strip(), "smiles", unique_key
    else:
        return "", "name", unique_key  # Return empty string and indicate it's from a name input

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

# Function to detect functional groups
def detect_functional_groups(mol):
    if not mol:
        return []
    
    functional_groups = []
    
    # Define SMARTS patterns for common functional groups
    patterns = {
        'Alcohol': '[OX2H]',
        'Aldehyde': '[CX3H1](=O)[#6]',
        'Ketone': '[#6][CX3](=O)[#6]',
        'Carboxylic acid': '[CX3](=O)[OX2H1]',
        'Ester': '[#6][CX3](=O)[OX2][#6]',
        'Amide': '[NX3][CX3](=[OX1])[#6]',
        'Amine (primary)': '[NX3;H2][#6]',
        'Amine (secondary)': '[NX3;H1]([#6])[#6]',
        'Amine (tertiary)': '[NX3]([#6])([#6])[#6]',
        'Alkene': '[C]=[C]',
        'Alkyne': '[C]#[C]',
        'Aromatic': 'c1ccccc1',
        'Nitro': '[N+](=O)[O-]',
        'Nitrile': '[C]#[N]',
        'Halogen (F)': '[F]',
        'Halogen (Cl)': '[Cl]',
        'Halogen (Br)': '[Br]',
        'Halogen (I)': '[I]',
        'Thiol': '[SX2H]',
        'Ether': '[OD2]([#6])[#6]',
        'Sulfonic acid': '[SX4](=O)(=O)[OX2H]',
        'Phosphoric acid': '[PX4](=O)([OX2H])([OX2H])[OX2H]'
    }
    
    for name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            functional_groups.append(name)
    
    return functional_groups

# Function to predict reaction type based on reactants and products
def predict_reaction_type(reactant_mols, product_mols):
    if not reactant_mols or not product_mols:
        return "Unknown"
    
    # Extract functional groups
    reactant_groups = []
    for mol in reactant_mols:
        reactant_groups.extend(detect_functional_groups(mol))
    
    product_groups = []
    for mol in product_mols:
        product_groups.extend(detect_functional_groups(mol))
    
    # Simple reaction type prediction based on functional groups
    if 'Alcohol' in reactant_groups and 'Carboxylic acid' in reactant_groups and 'Ester' in product_groups:
        return "Esterification"
    elif 'Ester' in reactant_groups and 'Alcohol' in product_groups:
        return "Hydrolysis (Ester)"
    elif 'Ketone' in reactant_groups and 'Alcohol' in product_groups:
        return "Reduction (Ketone to Alcohol)"
    elif 'Alcohol' in reactant_groups and 'Ketone' in product_groups:
        return "Oxidation (Alcohol to Ketone)"
    elif 'Carboxylic acid' in reactant_groups and 'Amide' in product_groups:
        return "Amide Formation"
    elif 'Alkene' in reactant_groups and 'Halogen' in reactant_groups:
        return "Halogenation"
    elif 'Aromatic' in reactant_groups and 'Nitro' in product_groups:
        return "Nitration"
    elif 'Amine (primary)' in reactant_groups and 'Amide' in product_groups:
        return "Amidation"
    elif 'Aldehyde' in reactant_groups and 'Alcohol' in product_groups:
        return "Reduction (Aldehyde to Alcohol)"
    elif 'Alkene' in reactant_groups and any('Halogen' in fg for fg in product_groups):
        return "Addition (Halogen to Alkene)"
    elif 'Alkyne' in reactant_groups and 'Alkene' in product_groups:
        return "Reduction (Alkyne to Alkene)"
    else:
        return "Complex/Unknown Reaction"

# Function to suggest solvents based on reaction type and functional groups
def suggest_solvents(reaction_type, reactant_groups, product_groups):
    all_groups = reactant_groups + product_groups
    
    # Default solvents (safer options)
    default_solvents = ["Water", "Ethanol", "Isopropanol", "Ethyl acetate"]
    
    # Specific recommendations based on reaction type
    solvents_by_reaction = {
        "Esterification": ["Toluene", "Heptane", "2-Me-THF", "Ethyl acetate"],
        "Hydrolysis (Ester)": ["Water", "Ethanol/Water", "THF/Water"],
        "Reduction (Ketone to Alcohol)": ["Methanol", "Ethanol", "Isopropanol", "THF"],
        "Oxidation (Alcohol to Ketone)": ["Acetone", "Ethyl acetate", "Dichloromethane"],
        "Amide Formation": ["DMF", "Acetonitrile", "THF", "2-Me-THF"],
        "Halogenation": ["Dichloromethane", "Ethyl acetate", "Acetone"],
        "Nitration": ["Water", "Acetic acid"],
        "Amidation": ["DMF", "THF", "Acetonitrile", "2-Me-THF"],
        "Reduction (Aldehyde to Alcohol)": ["Methanol", "Ethanol", "THF"],
        "Addition (Halogen to Alkene)": ["Dichloromethane", "Chloroform", "Carbon tetrachloride"],
        "Reduction (Alkyne to Alkene)": ["Methanol", "Ethanol", "THF"]
    }
    
    # Polarity-based recommendations
    if any(fg in all_groups for fg in ['Carboxylic acid', 'Alcohol', 'Amine', 'Amide']):
        polar_solvents = ["Water", "Methanol", "Ethanol", "Isopropanol", "Acetone", "Acetonitrile", "DMSO", "DMF"]
        return polar_solvents
    elif any(fg in all_groups for fg in ['Ester', 'Ketone', 'Aldehyde']):
        mid_solvents = ["Ethyl acetate", "THF", "2-Me-THF", "Acetone", "Dichloromethane"]
        return mid_solvents
    elif any(fg in all_groups for fg in ['Alkene', 'Alkyne', 'Aromatic', 'Halogen']):
        nonpolar_solvents = ["Toluene", "Heptane", "Cyclohexane", "Diethyl ether"]
        return nonpolar_solvents
    
    # If reaction type is known, use specific recommendations
    if reaction_type in solvents_by_reaction:
        return solvents_by_reaction[reaction_type]
    
    # Fall back to default solvents
    return default_solvents

# Function to suggest catalysts and reagents based on reaction type
def suggest_reagents(reaction_type, reactant_groups, product_groups):
    # Default reagents
    default_reagents = ["No specific reagent recommendation"]
    
    # Specific recommendations based on reaction type
    reagents_by_reaction = {
        "Esterification": ["Sulfuric acid", "p-Toluenesulfonic acid", "Amberlyst"],
        "Hydrolysis (Ester)": ["Sodium hydroxide", "Potassium hydroxide", "Hydrochloric acid"],
        "Reduction (Ketone to Alcohol)": ["Sodium borohydride", "Lithium aluminum hydride", "Hydrogen/Raney Ni"],
        "Oxidation (Alcohol to Ketone)": ["Pyridinium chlorochromate (PCC)", "Jones reagent", "TEMPO/NaOCl"],
        "Amide Formation": ["1-Ethyl-3-(3-dimethylaminopropyl)carbodiimide (EDC)", "Dicyclohexylcarbodiimide (DCC)"],
        "Halogenation": ["N-Bromosuccinimide (NBS)", "Bromine", "Chlorine"],
        "Nitration": ["Nitric acid", "Sulfuric acid"],
        "Amidation": ["EDC/HOBt", "DCC/NHS", "T3P"],
        "Reduction (Aldehyde to Alcohol)": ["Sodium borohydride", "Lithium aluminum hydride"],
        "Addition (Halogen to Alkene)": ["Bromine", "Chlorine", "Iodine"],
        "Reduction (Alkyne to Alkene)": ["Hydrogen/Lindlar catalyst", "Sodium in liquid ammonia"]
    }
    
    # If reaction type is known, use specific recommendations
    if reaction_type in reagents_by_reaction:
        return reagents_by_reaction[reaction_type]
    
    # Fallback based on functional groups
    if 'Carboxylic acid' in reactant_groups and 'Ester' in product_groups:
        return ["Sulfuric acid", "p-Toluenesulfonic acid"]
    elif 'Alcohol' in reactant_groups and 'Ketone' in product_groups:
        return ["PCC", "TEMPO/NaOCl", "Jones reagent"]
    elif 'Ketone' in reactant_groups and 'Alcohol' in product_groups:
        return ["Sodium borohydride", "Lithium aluminum hydride"]
    
    # Fall back to default reagents
    return default_reagents

# Function to suggest greener alternatives
def suggest_green_alternatives(solvents, reagents):
    green_solvents = {
        "Dichloromethane": ["2-Methyltetrahydrofuran", "Cyclopentyl methyl ether", "Ethyl acetate"],
        "Chloroform": ["2-Methyltetrahydrofuran", "Ethyl acetate", "Propylene carbonate"],
        "Carbon tetrachloride": ["D-Limonene", "2-Methyltetrahydrofuran"],
        "DMF": ["Ethyl lactate", "Propylene carbonate", "Cyrene"],
        "DMSO": ["Ethyl lactate", "Propylene carbonate"],
        "Benzene": ["Toluene", "Cyclohexane", "Heptane"],
        "Hexane": ["Heptane", "Cyclopentyl methyl ether"]
    }
    
    green_reagents = {
        "Jones reagent": ["TEMPO/NaOCl", "Oxone"],
        "PCC": ["TEMPO/NaOCl", "IBX in Ethyl lactate"],
        "Chromium oxide": ["TEMPO/NaOCl", "IBX"],
        "n-Butyllithium": ["Turbo Grignard", "Activated Mg"],
        "Trifluoroacetic acid": ["Acetic acid", "Formic acid"],
        "Sodium azide": ["1H-Tetrazole derivatives"] 
    }
    
    green_solvent_alternatives = []
    for solvent in solvents:
        if solvent in green_solvents:
            green_solvent_alternatives.append(f"{solvent} → {', '.join(green_solvents[solvent])}")
    
    green_reagent_alternatives = []
    for reagent in reagents:
        if reagent in green_reagents:
            green_reagent_alternatives.append(f"{reagent} → {', '.join(green_reagents[reagent])}")
    
    return green_solvent_alternatives, green_reagent_alternatives

# Function to calculate molecular properties
def calculate_properties(mol):
    if not mol:
        return {}
    
    properties = {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2)
    }
    
    return properties

# Container for inputs
with st.container():
    st.header("Reactants")
    reactants_data = []
    
    # Dynamic reactants
    for i in range(1, st.session_state.reactants_count + 1):
        col1, col2 = st.columns([10, 1])
        with col1:
            reactant_result = get_smiles_input("Reactant", "reactant", i)
            reactants_data.append(reactant_result)
        with col2:
            if i == st.session_state.reactants_count and i > 1:
                st.button("➖", key=f"remove_reactant_{i}", on_click=remove_component, args=('reactant',))
    
    # Button to add new reactant
    st.button("Add Reactant ➕", on_click=add_component, args=('reactant',), key="add_reactant")

with st.container():
    st.header("Products")
    products_data = []
    
    # Dynamic products
    for i in range(1, st.session_state.products_count + 1):
        col1, col2 = st.columns([10, 1])
        with col1:
            product_result = get_smiles_input("Product", "product", i)
            products_data.append(product_result)
        with col2:
            if i == st.session_state.products_count and i > 1:
                st.button("➖", key=f"remove_product_{i}", on_click=remove_component, args=('product',))
    
    # Button to add new product
    st.button("Add Product ➕", on_click=add_component, args=('product',), key="add_product")

with st.container():
    st.header("Agents (optional)")
    agents_data = []
    
    # Dynamic agents
    for i in range(1, st.session_state.agents_count + 1):
        col1, col2 = st.columns([10, 1])
        with col1:
            agent_result = get_smiles_input("Agent", "agent", i)
            agents_data.append(agent_result)
        with col2:
            if i == st.session_state.agents_count and i > 1:
                st.button("➖", key=f"remove_agent_{i}", on_click=remove_component, args=('agent',))
    
    # Button to add new agent
    st.button("Add Agent ➕", on_click=add_component, args=('agent',), key="add_agent")

# Submit and show results
if st.button("Submit", key="submit_btn"):
    with st.spinner("Processing chemical data and making predictions..."):
        # Process all components
        processed_reactants = []
        processed_products = []
        processed_agents = []
        
        # Process reactants
        for smiles, source, unique_key in reactants_data:
            if source == "name" and f"{unique_key}_chem_name" in st.session_state:
                name = st.session_state[f"{unique_key}_chem_name"]
                converted_smiles = name_to_smiles(name)
                if converted_smiles:
                    st.success(f"Successfully converted '{name}' to SMILES")
                    processed_reactants.append({
                        "smiles": converted_smiles,
                        "source": "name", 
                        "name": name,
                        "key": unique_key,
                        "mol": Chem.MolFromSmiles(converted_smiles)
                    })
                else:
                    st.error(f"Could not find SMILES for '{name}'")
            elif smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    processed_reactants.append({
                        "smiles": smiles, 
                        "source": source,
                        "key": unique_key,
                        "mol": mol
                    })
                else:
                    st.error(f"Invalid SMILES: {smiles}")
        
        # Process products
        for smiles, source, unique_key in products_data:
            if source == "name" and f"{unique_key}_chem_name" in st.session_state:
                name = st.session_state[f"{unique_key}_chem_name"]
                converted_smiles = name_to_smiles(name)
                if converted_smiles:
                    st.success(f"Successfully converted '{name}' to SMILES")
                    processed_products.append({
                        "smiles": converted_smiles,
                        "source": "name", 
                        "name": name,
                        "key": unique_key,
                        "mol": Chem.MolFromSmiles(converted_smiles)
                    })
                else:
                    st.error(f"Could not find SMILES for '{name}'")
            elif smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    processed_products.append({
                        "smiles": smiles, 
                        "source": source,
                        "key": unique_key,
                        "mol": mol
                    })
                else:
                    st.error(f"Invalid SMILES: {smiles}")
        
        # Process agents
        for smiles, source, unique_key in agents_data:
            if source == "name" and f"{unique_key}_chem_name" in st.session_state:
                name = st.session_state[f"{unique_key}_chem_name"]
                converted_smiles = name_to_smiles(name)
                if converted_smiles:
                    st.success(f"Successfully converted '{name}' to SMILES")
                    processed_agents.append({
                        "smiles": converted_smiles,
                        "source": "name", 
                        "name": name,
                        "key": unique_key,
                        "mol": Chem.MolFromSmiles(converted_smiles)
                    })
                else:
                    st.error(f"Could not find SMILES for '{name}'")
            elif smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    processed_agents.append({
                        "smiles": smiles, 
                        "source": source,
                        "key": unique_key,
                        "mol": mol
                    })
                else:
                    st.error(f"Invalid SMILES: {smiles}")
        
        # Short delay to allow the success/error messages to be displayed
        time.sleep(0.5)

    # Prepare reaction data
    reaction_data = {
        "reactants": [item["smiles"] for item in processed_reactants],
        "products": [item["smiles"] for item in processed_products],
        "agents": [item["smiles"] for item in processed_agents]
    }

    st.subheader("Parsed Reaction")
    st.json(reaction_data)

    # Show molecules if available
    if processed_reactants:
        st.subheader("Reactant Structures")
        for i, reactant in enumerate(processed_reactants, 1):
            with st.container():
                st.write(f"**Reactant {i}**")
                draw_mol(reactant["smiles"], f"Reactant {i}")
                if reactant.get("source") == "name":
                    st.caption(f"Generated from chemical name: {reactant['name']}")
                
                # Display molecular properties in an expander
                with st.expander("Molecular Properties"):
                    props = calculate_properties(reactant["mol"])
                    for prop, value in props.items():
                        st.write(f"**{prop}:** {value}")

    if processed_products:
        st.subheader("Product Structures")
        for i, product in enumerate(processed_products, 1):
            with st.container():
                st.write(f"**Product {i}**")
                draw_mol(product["smiles"], f"Product {i}")
                if product.get("source") == "name":
                    st.caption(f"Generated from chemical name: {product['name']}")
                
                # Display molecular properties in an expander
                with st.expander("Molecular Properties"):
                    props = calculate_properties(product["mol"])
                    for prop, value in props.items():
                        st.write(f"**{prop}:** {value}")

    if processed_agents:
        st.subheader("Agent Structures")
        for i, agent in enumerate(processed_agents, 1):
            with st.container():
                st.write(f"**Agent {i}**")
                draw_mol(agent["smiles"], f"Agent {i}")
                if agent.get("source") == "name":
                    st.caption(f"Generated from chemical name: {agent['name']}")
                
                # Display molecular properties in an expander
                with st.expander("Molecular Properties"):
                    props = calculate_properties(agent["mol"])
                    for prop, value in props.items():
                        st.write(f"**{prop}:** {value}")
    
    # Make predictions if we have both reactants and products
    if processed_reactants and processed_products:
        st.header("Reaction Analysis & Predictions")
        
        # Extract functional groups
        reactant_mols = [r["mol"] for r in processed_reactants]
        product_mols = [p["mol"] for p in processed_products]
        
        all_reactant_groups = []
        for mol in reactant_mols:
            groups = detect_functional_groups(mol)
            all_reactant_groups.extend(groups)
        
        all_product_groups = []
        for mol in product_mols:
            groups = detect_functional_groups(mol)
            all_product_groups.extend(groups)
        
        # Predict reaction type
        reaction_type = predict_reaction_type(reactant_mols, product_mols)
        
        st.subheader("Detected Reaction")
        st.write(f"**Reaction Type:** {reaction_type}")
        
        col1, col2 = st.columns(2)
        with col1:
            st.write("**Functional Groups in Reactants:**")
            if all_reactant_groups:
                for group in set(all_reactant_groups):
                    st.write(f"- {group}")
            else:
                st.write("- No specific functional groups detected")
        
        with col2:
            st.write("**Functional Groups in Products:**")
            if all_product_groups:
                for group in set(all_product_groups):
                    st.write(f"- {group}")
            else:
                st.write("- No specific functional groups detected")
        
        # Suggest solvents and reagents
        suggested_solvents = suggest_solvents(reaction_type, all_reactant_groups, all_product_groups)
        suggested_reagents = suggest_reagents(reaction_type, all_reactant_groups, all_product_groups)
        
        st.subheader("Recommended Reaction Conditions")
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**Suggested Solvents:**")
            for solvent in suggested_solvents:
                st.write(f"- {solvent}")
        
        with col2:
            st.write("**Suggested Reagents/Catalysts:**")
            for reagent in suggested_reagents:
                st.write(f"- {reagent}")
        
        # Green alternatives
        green_solvents, green_reagents = suggest_green_alternatives(suggested_solvents, suggested_reagents)
        
        with st.expander("Green Chemistry Alternatives"):
            if green_solvents:
                st.write("**Greener Solvent Alternatives:**")
                for alt in green_solvents:
                    st.write(f"- {alt}")
            else:
                st.write("- The suggested solvents are already considered green options")
            
            if green_reagents:
                st.write("**Greener Reagent Alternatives:**")
                for alt in green_reagents:
                    st.write(f"- {alt}")
            else:
                st.write("- The suggested reagents are already considered green options")
"""
Lookup MW, toxicity, GHS, flashpoint, solvent “greenness” via PubChem/CHEM21.
"""

import requests
import re
import os
import streamlit as st

@st.cache_data
def get_cid(smiles: str)->str:
    # Get the CID from SMILES
    # A conversion is made from SMILES to CID (Compound ID) because
    # PubChem is organised with unique identifiers (CID) for each known molecule.
    cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    cid_response = requests.get(cid_url)
   
    if cid_response.status_code != 200:
        return "Error when retrieving CID."
   
    # When a request is made to the PubChem API, the response is a JSON (JavaScript Object Notation) object.
    # so cid_response.json() retrieves the JSON response from the PubChem API.
    # Next, .get(‘IdentifierList’, {}) accesses the IdentifierList key in the JSON response, which contains the list of CIDs.
    # Finally, .get(‘CID’, []) is used to extract the list of CIDs.
   
    cids = cid_response.json().get("IdentifierList", {}).get("CID", [])
    if not cids:
        return "No CID found for this SMILES."
   
    # In PubChem, a molecule can have several CIDs associated with it,  
    # cids[0] selects the first identifier in the list
   
    cid = cids[0]
    return cid

@st.cache_data
def get_flash_point_from_smiles(smiles: str):

    cid=get_cid(smiles)

    # Step 1 : Obtenir les propriétés expérimentales (dont le flash point)
    prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Experimental+Properties"
    prop_response = requests.get(prop_url)

    if prop_response.status_code != 200:
        return "Error retrieving properties"

    # Étape 2 : Chercher "Flash Point" dans les propriétés
    data = prop_response.json()
    sections = data.get("Record", {}).get("Section", [])

    for section in sections:
        if section.get("TOCHeading")=="Chemical and Physical Properties":
            for subsection in section.get("Section", []):
                if subsection.get("TOCHeading") == "Experimental Properties":
                    for sub_subsection in subsection.get("Section", []):
                        if sub_subsection.get("TOCHeading") == "Flash Point":
                            infos=sub_subsection.get("Information", [])
                            info=infos[0]
                            # info.get("Value", {})
                            # info est un dictionnaire dans lequel chaque élément de la réponse JSON est stocké
                            # "Value" est une clé qui contient généralement des informations sur une propriété spécifique de la molécule
                            # comme le flashpoint. Si "Value" n'est pas présent un dictionnaire vide ({}) est renvoyé pour éviter les erreurs.
                           
                            # .get("StringWithMarkup", [{}])
                            # "StringWithMarkup" est un champ spécifique dans la structure JSON retournée par PubChem, et il contient généralement
                            # une liste d'objets qui possèdent un champ String avec la valeur qu'on recherche.
                            # "StringWithMarkup" est généralement une liste d'objets, et on doit accéder à cet objet avec un indice pour ça
                            # qu'il y a [{}].
                           
                            # [0].get("String")
                            # Comme "StringWithMarkup" est une liste d'objets, on veut son premier élément.
                            # À ce moment, cette partie est un dictionnaire et on veut obtenir la valeur de la clé "String",
                            # qui contient le texte brut de la propriété.
                            value = info.get("Value", {}).get("StringWithMarkup", [{}])[0].get("String")
                           
                            # Expression régulière pour extraire la température et l'unité (F ou C)
                            # La fonction re.match() sert à extraire la température et l'unité à partir de la chaîne de texte.
                            # L'expression r"([+-]?\d+\.?\d*)\s*°?\s*(F|C)" est conçue pour capturer
                            # un nombre entier ou décimal (([+-]?\d+\.?\d*)) et le symbol de l'unité.
                            match = re.match(r"([+-]?\d+\.?\d*)\s*°?\s*(F|C)", value)
                            if match:
                             # Extraire la valeur de température et l'unité
                                temp_value = float(match.group(1))
                                unit = match.group(2).upper()
       
                                if unit == "F":
                                # Si la température est en Fahrenheit, on la convertit en Celsius
                                    celsius = (temp_value - 32) * 5.0 / 9.0
                                    return round(celsius, 2)
       
                                elif unit == "C":
                                # Si la température est déjà en Celsius, on la retourne directement
                                    return round(temp_value, 2)
   
                            return "Temperature in invalid format"
                                   
    return "Flash point not found for this molecule."

@st.cache_data
def get_ghs_data(smiles:str):
    cid=get_cid(smiles)

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    response = requests.get(url)
    if response.status_code != 200:
        return "Error retrieving properties"

    data = response.json()
    ghs_info = []
   
    # L'information GHS est dans la section "GHS Classification"
    for section in data['Record']['Section']:
        if section['TOCHeading'] == 'Safety and Hazards':
            for sub_section in section['Section']:
                if sub_section['TOCHeading'] == 'Hazards Identification':
                    for sub_subsection in sub_section["Section"]:
                        if sub_subsection['TOCHeading'] == 'GHS Classification':
                            for info in sub_subsection['Information']:
                                if 'Value' in info and 'StringWithMarkup' in info['Value']:
                                    for item in info['Value']['StringWithMarkup']:
                                        text = item['String']
                                        ghs_info.append(text)
                                        h_codes = []
                                        blank_line_count = 0
                                        for line in ghs_info:
                                            if line.strip() == '':  
                                                blank_line_count += 1
                                                if blank_line_count == 2:
                                                    break
                                                # Find H-codes in
                                            else:
                                                found = re.findall(r'H\d{3}', line)
                                                h_codes.extend(found)
                                 
    return h_codes

@st.cache_data
def hazard_statements(smiles: str):
    hazard_pictograms=[]
    explosives=["H200","H201","H202","H203","H204","H209","H210","H211","H240","H241","H280","H284"]
    flammables=["H205","H206","H207","H208","H220","H221","H222","H223","H224","H225",
                "H226","H228","H229","H230","H231","H232","H241","H242","H250","H251",
                "H252","H260","H261","H282","H283"]
    oxidizers= ["H270","H271","H272"]
    corrosives=["H290","H314","H318"]
    acute_toxicity=["H300","H301","H310","H311","H330","H331"]
    health_hazard=["H304","H305","H334","H340","H341","H350","H351","H360","H361","H370","H371","H372","H373"]
    environment_hazard=["H400","H401","H402","H410","H411","H412","H413","H420"]
    h_codes=get_ghs_data(smiles)
    #
    for h_code in h_codes:
        if h_code in explosives:
            hazard_pictograms.append("Explosive")
        elif h_code in flammables:
            hazard_pictograms.append("Flammable")
        elif h_code in oxidizers:
            hazard_pictograms.append("Oxidizing")
        elif h_code in corrosives:
            hazard_pictograms.append("Corrosive")
        elif h_code in acute_toxicity:
            hazard_pictograms.append("Toxic")
        elif h_code in health_hazard:
            hazard_pictograms.append("Health Hazard")
        elif h_code in environment_hazard:
            hazard_pictograms.append("Environment Hazard")
    #
    unique_hazard_pictograms=list(set(hazard_pictograms))
    return unique_hazard_pictograms

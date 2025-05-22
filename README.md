![Project Logo](assets/banner2.jpeg)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
MCGA
</h1>

<br>

# MCGA

A Python toolkit to score the greenness of chemical reactions by the 12 Principles of Green Chemistry.

---

## 📖 Table of Contents

- [🧪 About](#-about)
- [🚀 Features](#-features)
- [👩‍💻 Installation](#-installation)
- [🚀 Quick Start](#-quick-start)
- [🧪 Testing](#-testing)
- [♻️ 12 Principles Check-List](#-12-principles-check-list)
- [📜 License](#-license)
- [👥 Contributors](#-contributors)

---
## 🧪 About

MCGA (**M**ake **C**hemistry **G**reen **A**gain) is a toolkit for evaluating chemical reactions according to the 12 Principles of Green Chemistry.  
It automates the scoring, hazard identification, and reaction visualization, providing both a Streamlit web interface and a Python package.

---

## 🚀 Features

- **12-Principles Checklist:** Quantifies how “green” a reaction is.
- **Automated Metrics:** Calculates Atom Economy, E-Factor, and other green chemistry indicators.
- **Hazard Identification:** GHS hazard code extraction, flash-point prediction, and pictograms.
- **Solvent & Catalyst Prediction:** via Gemini API (optional)
- **Streamlit Web App:** Interactive interface for entering reactions and visualizing results.
- **Flexible Input:** Accepts chemical names, SMILES, or drawn structures.

---

## 👩‍💻 Installation

**Prerequisites:**  
- Python 3.11+  
- [Conda](https://docs.conda.io/en/latest/) recommended  
- [RDKit](https://www.rdkit.org/) (will install via conda)  

**To create a new environment and install MCGA:**
```bash
conda create -n mcga python=3.11
conda activate mcga
conda install -c conda-forge rdkit
pip install -e .
```

---

## 🚀 Quick Start

First, open a terminal and make sure you are in the correct folder:

```
(conda_env) $ cd yourpathto/MCGA/src/mcga/ # Replace with your actual path; make sure you're in 'mcga/'
(conda_env) $ streamlit run app.py
```
After running this command, Streamlit will start the app and open it in your web browser. If your browser doesn’t open automatically, copy and paste the local URL provided in the terminal into your browser’s address bar.

**What you will see**
On the screen, you will see three main modules (Reactants, Product, Agents), where you can enter one or more compounds using their chemical name, SMILES string, or by drawing the structure. Each module allows you to add more components by clicking the Add Reactant, Add Product, or Add Agent button.

Once all your components are entered, click Submit reaction to evaluate your reaction using the Green Chemistry criteria.

---

## 🧪 Testing
To run all tests and check coverage:
```
(conda_env) $ pip install tox
(conda_env) $ tox
```
Test files are in the tests/ folder and cover all core functionality.


## ♻️ 12 Principles Check-List

| Principle                                          | “Green” if…                                                                                                   |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| 1. Prevention of Waste                             | E-factor < 25 or PMI < 100                                                                                     |
| 2. Atom Economy                                    | % Atom Economy ≥ 90 %                                                                                         |
| 3. Less Hazardous Syntheses                        | No GHS Acute Tox. Cat. 1–2 or CMR Cat. 1–2 among reagents/by-products                                          |
| 4. Designing Safer Chemicals                       | Final product LD₅₀ (oral, rat) > 2 000 mg/kg; no ≥ Cat 1 aquatic toxicity                                       |
| 5. Safer Solvents & Auxiliaries                    | All solvents from CHEM21 “recommended”; auxiliaries ≤ 10 % w/w                                                 |
| 6. Design for Energy Efficiency                    | Reaction temperature ≤ 50 °C & ambient pressure only                                                          |
| 7. Use of Renewable Feedstocks                     | ≥ 50 % of total carbon atoms from bio-based feedstocks                                                        |
| 8. Reduce Derivatives                              | ≤ 1 protection/deprotection step                                                                               |
| 9. Catalysis                                       | Uses catalyst loading ≤ 10 mol %                                                                               |
| 10. Design for Degradation                         | Predicted environmental half-life < 60 days; no PBT (persistent/bioaccumulative/toxic) flags                  |
| 11. Real-Time Analysis for Pollution Prevention    | At least one in-line monitor (FTIR, GC, HPLC…)                                                                  |
| 12. Inherently Safer Chemistry (Accident Prevention) | All reagents flash-point ≥ 60 °C; no peroxides or explosophoric groups                                         |

## 📜 License

This project is licensed under the MIT License

## 👥 Contributors

Team Members and Main Roles

Jiho Keum (@jihokeum) — Streamlit UI, Harazdous by-product module, fire/explosion module, README file, 

Alexia Dade (@alexiadade) —

Bilel Bouzouaid (@BilelBouzouaid) —

Ylann Willemin (@Ylann-Willemin) -

All members contributed to coding, and testing.




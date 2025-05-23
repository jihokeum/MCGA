![Project Logo](assets/banner2.jpeg)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
MCGA
</h1>

<br>

A Python toolkit to score the greenness of chemical reactions by the 12 Principles of Green Chemistry.


## 📖 Table of Contents

- [🧪 About](#-about)
- [👩‍💻 Installation](#-installation)
- [🚀 Quick Start](#-quick-start)
- [🧪 Testing](#-testing)
- [♻️ 12 Principles Check-List](#-12-principles-check-list)
- [📜 License](#-license)
- [👥 Contributors](#-contributors)


## 🧪 About

MCGA (**M**ake **C**hemistry **G**reen **A**gain) is a toolkit for evaluating chemical reactions according to the 12 Principles of Green Chemistry.  
It automates the scoring, hazard identification, and reaction visualization, providing both a Streamlit web interface and a Python package.
It automates reaction scoring, hazard identification, and visualization through a user-friendly Streamlit web interface and Python package.

In this project, we specifically focused on:

* Accident prevention (fire and explosion risk)

* Health and environmental hazards (toxicity of by-products)

* Waste reduction (e-factor)

* Optimization of desired product (atom economy)

The thresholds used for “green” (such as minimum flash point, maximum E-factor, etc.) were determined from the scientific literature and group discussions. These cutoffs are indicative and should not be considered strict or universal standards. Other Green Chemistry principles are not evaluated in this program.


## 👩‍💻 Installation

**Prerequisites:**  
- Python 3.11.x  
- [Conda](https://docs.conda.io/en/latest/) recommended  
- [RDKit](https://www.rdkit.org/) (see below)

**To create a new environment and install MCGA:**
```bash
git clone https://github.com/jihokeum/MCGA
cd yourpathto/MCGA/ # Replace with your actual path
conda create -n mcga python=3.11
conda activate mcga
conda install -c conda-forge rdkit
pip install -e .
```


## 🚀 Quick Start

First, open a terminal and make sure you are in the correct folder:

```
(conda_env) $ cd yourpathto/MCGA/src/mcga/ # Replace with your actual path
(conda_env) $ streamlit run app.py
```
After running this command, Streamlit will start the app and open it in your web browser. If your browser doesn’t open automatically, copy and paste the local URL provided in the terminal into your browser’s address bar.

#### ➡️ Example Interface 1

<img src="assets/interface1.png" alt="Interface 1" width="500"/>

On the screen, you will see three main modules (Reactants, Product, Agents), where you can enter one or more compounds using their chemical name, SMILES string, or by drawing the structure. Each module allows you to add more components. Once all your components are entered, click **Submit reaction** to evaluate your reaction using the Green Chemistry criteria.

#### ➡️ Example Interface 2
<img src="assets/interface2.png" alt="Interface 2" width="500"/>
<img src="assets/interface3.png" alt="Interface 3" width="500"/>
<img src="assets/interface4.png" alt="Interface 4" width="500"/>


## 🧪 Testing
To run all tests and check coverage:
```
(conda_env) $ pip install pytest
(conda_env) $ pytest 
```
Test files are in the tests/ folder and cover all core functionality.


## ♻️ 12 Principles Check-List

MCGA evaluates the following Green Chemistry Principles:

| Principle                                          | “Green” if…                                                                                                   |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| 1. Prevention of Waste                             | E-factor < 5 (excellent), < 25 (moderate)                                                                                     |
| 2. Atom Economy                                    | % Atom Economy ≥ 75 % (excellent),  ≥ 50 % (moderate)                                                                                         |
| 3. Less Hazardous Syntheses                        | No by-product with GHS acute toxicity or CMR code (excellent), 1-3 GHS acute toxicity code (moderate)                                         |
| 12. Inherently Safer Chemistry (Accident Prevention) | No peroxides or explosophoric groups, and all reagents flash-point ≥ 60 °C (excellent), ≥ 20 °C (moderate)                                      |


## 📜 License

This project is licensed under the MIT License


## 👥 Contributors

Team Members and Main Roles

Jiho Keum (@jihokeum) — Streamlit UI, Full reaction scheme module, Hazardous by-product module, fire/explosion module, Compound details module, README file, test files

Alexia Dade (@alexiadade) — chemistry green list interface, traffic light grading, AE/E-factor module, balancing equations (balancing_equations.py), Streamlit UI

Bilel Bouzouaid (@BilelBouzouaid) — condition prediction via Gemini, Streamlit UI, AE/E-factor module,Compound details module

Ylann Willemin (@Ylann-Willemin) - getting physical properties via Pubchem (lookup.py)

All members contributed to coding, writing the notebook, searching and finding useful packages.




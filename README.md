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
conda create -n mcga python=3.11
conda activate mcga
conda install -c conda-forge rdkit
pip install -e .
```


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


## 🧪 Testing
To run all tests and check coverage:
```
(conda_env) $ pip install tox
(conda_env) $ tox
```
Test files are in the tests/ folder and cover all core functionality.


## ♻️ 12 Principles Check-List

Currently, MCGA automatically evaluates the following Green Chemistry Principles:

| Principle                                          | “Green” if…                                                                                                   |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| 1. Prevention of Waste                             | E-factor < 5 (excellent), < 25(moderate)                                                                                     |
| 2. Atom Economy                                    | % Atom Economy ≥ 75 % (excellent),  ≥ 50 % (moderate)                                                                                         |
| 3. Less Hazardous Syntheses                        | No by-product with more then 4 GHS acute toxicity codes or at least one CMR code                                          |
| 12. Inherently Safer Chemistry (Accident Prevention) | All reagents flash-point ≥ 60 °C; no peroxides or explosophoric groups (excellent)                                        |


## 📜 License

This project is licensed under the MIT License


## 👥 Contributors

Team Members and Main Roles

Jiho Keum (@jihokeum) — Streamlit UI, Harazdous by-product module, fire/explosion module, README file, 

Alexia Dade (@alexiadade) —

Bilel Bouzouaid (@BilelBouzouaid) —

Ylann Willemin (@Ylann-Willemin) -

All members contributed to coding, and testing.




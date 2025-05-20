![Project Logo](assets/banner2.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
MCGA
</h1>

<br>

## 🔥 Usage

```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

## 👩‍💻 Installation

To initialize repo:
Create a new environment, you may also give the environment a different name. 

```
conda create -n mcga python=3.11
conda activate mcga
pip install -e .
```

To run the program:
```
cd MCGA/src/mcga
streamlit run app.py
```

# MCGA

A Python toolkit to score the greenness of chemical reactions by the 12 Principles of Green Chemistry.

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

## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:jihokeum/mcga`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:jihokeum/mcga.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(mcga) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```




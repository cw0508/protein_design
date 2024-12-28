# Protein Design Projects

This repository contains solutions and analyses for assignments related to **Protein Design** as part of the CSCI-GA.3033-111 course at New York University. Each assignment explores different tools and methodologies for computational protein design, prediction, and analysis.

---

## Table of Contents

1. [Assignment Summaries](#assignment-summaries)
2. [Technologies and Tools](#technologies-and-tools)
3. [Usage Instructions](#usage-instructions)
4. [Results Overview](#results-overview)
5. [Acknowledgments](#acknowledgments)

---

## Assignment Summaries

### **HW1: Protein Data Exploration**

- Explored the structure and function of Trypsin using UniProt and PDB.
- Performed sequence and structure similarity searches.
- Visualized and aligned protein structures using PyMol.
- Implemented sequence and structural distance functions in Python.
- Benchmarked protein sets for sequence and structure similarity.

### **HW2: Protein Structure Prediction**

- Predicted structures for proteins using AlphaFold2, AlphaFold3, and RosettaFold.
- Generated pLDDT and PAE charts to evaluate model performance.
- Aligned predicted structures with experimental ones and compared key differences.
- Investigated the prediction of ECORI with and without DNA binding sites.

### **HW3: Binder Design with BindCraft and Chai-1**

- Used BindCraft to design binders for the BHRF1 protein targeting specific hotspots.
- Validated binders using Chai-1 and compared predicted structures.
- Evaluated in-silico success metrics such as RMSD, Binding Energy, SASA, and pLDDT.
- Visualized passing and non-passing binders to analyze performance differences.

### **HW4: Protein Design with RFDiffusion**

- Generated protein designs using RFDiffusion, including symmetric homo-trimers and motif scaffolds.
- Visualized alpha-helix fractions and compared scaffolds with partial diffusion outputs.
- Designed and evaluated binders for IL-7Rα, optimizing hotspots and structural stability.
- Analyzed protein-ligand complexes using PyRosetta for energy calculations.

---

## Technologies and Tools

- **UniProt** and **PDB**: Protein data exploration.
- **PyMol**: Visualization and alignment of protein structures.
- **AlphaFold2/3**, **RosettaFold**, and **RFDiffusion**: Protein structure prediction and design.
- **BindCraft** and **Chai-1**: Binder design and validation.
- **PyRosetta**: Energy calculations for protein-protein and protein-ligand complexes.
- **Python Libraries**:
  - `Bio.PDB`: Sequence extraction and structural analysis.
  - `matplotlib` and `numpy`: Visualization and data analysis.
  - Custom implementations for sequence and structure distance metrics.

---

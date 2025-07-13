# Data Cleaning Script for USPTO-all-remapped Used in ChemEnzyRetroPlanner Single-Step Model

## Data Preparation

The files in `data/USPTO_remapped` are downloaded from [IBM RXNMapperData](https://ibm.ent.box.com/v/RXNMapperData/folder/112951098080).

## Packages
- rdchiral
- pandas 
- rdkit

## Pipeline Scripts

This repository provides a pipeline for cleaning, processing, and extracting reaction templates from USPTO reaction data. The main scripts are:

### 1. `1.clean_remaped_uspto.py`
- **Function:** Cleans the raw mapped USPTO reaction files by removing unmapped reactants and ensuring atom mapping consistency between reactants and products.
- **Input:** Files in `data/USPTO_remapped` (downloaded from the above link).
- **Output:** Cleaned reaction data saved to `data/Cleaned_data/all_clean_remapped_USPTO.csv`.

### 2. `2.extract_retrotemplates.py`
- **Function:** Extracts retrosynthetic templates from the cleaned reaction data.
- **Input:** `data/Cleaned_data/all_clean_remapped_USPTO.csv`
- **Output:** Extracted templates saved to `data/Cleaned_data/USPTO_remapped_rxn_templates.csv`.

### 3. `3.remove_same_between_react_and_prod.py`
- **Function:** Removes reactions where the same molecule appears in both reactants and products, further deduplicating the dataset.
- **Input:** `data/Cleaned_data/USPTO_remapped_rxn_templates.csv`
- **Output:** Final deduplicated templates saved to `data/Cleaned_data/USPTO_remapped_remove_same_rxn_templates.csv`

## Directory Structure

```
project_root/
├── data/
│   ├── USPTO_remapped/         # Downloaded raw mapped USPTO files
│   └── Cleaned_data/           # Outputs from each pipeline step
├── 1.clean_remaped_uspto.py
├── 2.extract_retrotemplates.py
├── 3.remove_same_between_react_and_prod.py
└── README.md
```

## Reference
- [IBM RXNMapperData Download Link](https://ibm.ent.box.com/v/RXNMapperData/folder/112951098080) 
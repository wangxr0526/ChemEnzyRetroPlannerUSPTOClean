import re

import numpy as np
import pandas as pd
from tqdm import tqdm

from utils import canonicalize_smiles

def remove_same_mol(rxn):
    pt = re.compile(r':(\d+)]')
    reacts, prods = rxn.split('>>')
    is_remove = False
    if '.' in prods:
        reacts_list = reacts.split('.')
        prods_list = np.asanyarray(prods.split('.'))
        new_reacts_list = []
        for react in reacts_list:
            if (react == prods_list).sum() > 0:
                is_remove = True
                prods_list = prods_list[(prods_list != react)]
            else:
                new_reacts_list.append(react)
        if is_remove:
            prods_list = prods_list.tolist()
            new_reacts = canonicalize_smiles('.'.join(new_reacts_list), clear_map=False)
            map_new_reacts = sorted(re.findall(pt, new_reacts))
            new_prods = canonicalize_smiles('.'.join(prods_list), clear_map=False)
            map_new_prods = sorted(re.findall(pt, new_prods))
            if map_new_prods == map_new_reacts:
                rxn = '{}>>{}'.format(new_reacts, new_prods)
                if rxn not in ['>>', '.>>.']:
                    return rxn, is_remove
                else:
                    return None, is_remove
            else:
                print(f'{rxn} bad map after remove same mol!!!')
                return rxn, False
        else:
            return rxn, is_remove
    else:
        return rxn, is_remove
    



if __name__ == '__main__':
    extract_df = pd.read_csv('data/Cleaned_data/USPTO_remapped_rxn_templates.csv')

    patent_ids = extract_df['source'].tolist()
    droped_unmapped_rxns = extract_df['droped_unmapped_rxn'].tolist()
    retro_templates = extract_df['retro_template'].tolist()
    clean_map_rxns = extract_df['clean_map_rxn'].tolist()

    empty_after_remove = 0
    unique_rxn_dict = {}
    for patent_id, rxn, template, clean_map_rxn in tqdm(
            zip(patent_ids, droped_unmapped_rxns, retro_templates, clean_map_rxns),
            total=len(droped_unmapped_rxns), desc='去除产物与反应物相同分子'):
        rxn, _ = remove_same_mol(rxn)
        clean_map_rxn, _ = remove_same_mol(clean_map_rxn)
        if rxn and clean_map_rxn is not None:
            if clean_map_rxn not in unique_rxn_dict:
                unique_rxn_dict[clean_map_rxn] = (patent_id, rxn, template)
        else:
            empty_after_remove += 1
            print(f'{patent_id} 去除后反应为空: {rxn}')
    result_patent_ids = []
    result_rxns = []
    result_templates = []
    result_clean_map_rxns = []
    for clean_map_rxn, (patent_id, rxn, template) in list(unique_rxn_dict.items()):
        result_patent_ids.append(patent_id)
        result_rxns.append(rxn)
        result_templates.append(template)
        result_clean_map_rxns.append(clean_map_rxn)
    print(f'去除后为空的反应数: {empty_after_remove}，成功保留: {len(result_patent_ids)} 条')
    result_df = pd.DataFrame({})
    result_df['source'] = result_patent_ids
    result_df['droped_unmapped_rxn'] = result_rxns
    result_df['retro_template'] = result_templates
    result_df['clean_map_rxn'] = result_clean_map_rxns
    
        
    uspto_190_df = pd.read_csv('data/USPTO-multistep-190/uspto_190_route_target_mol.csv')
    target_products = set(uspto_190_df['smiles'])  

    # 去除产物在 uspto_190 中的反应
    result_df = result_df[result_df['clean_map_rxn'].apply(
        lambda rxn: rxn.split('>>')[-1] not in target_products
    )].reset_index(drop=True)
    
    result_df.to_csv('data/Cleaned_data/USPTO_remapped_remove_same_rxn_templates.csv', index=False)

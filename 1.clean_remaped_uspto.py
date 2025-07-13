import os
import re
import multiprocessing
import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem
from tqdm import tqdm

from utils import canonicalize_smiles

debug = False

def drop_unmapped_reactants(rxn):
    pt = re.compile(r':(\d+)]')
    reacts, _, prod = rxn.split('>')
    react_list = reacts.split('.')
    new_react_list = []
    for react in react_list:
        if re.findall(pt, react):
            new_react_list.append(react)
    new_reacts = '.'.join(new_react_list)
    react_maps = sorted(re.findall(pt, new_reacts))
    prod_maps = sorted(re.findall(pt, prod))
    clean_map_rxn = '{}>>{}'.format(canonicalize_smiles(new_reacts, clear_map=True),
                                    canonicalize_smiles(prod, clear_map=True))
    drop_unmapped_rxn = '{}>>{}'.format(canonicalize_smiles(new_reacts, clear_map=False),
                                        canonicalize_smiles(prod, clear_map=False))
    if react_maps != prod_maps:
        bad_map_flag = True
    else:
        bad_map_flag = False
    return drop_unmapped_rxn, clean_map_rxn, bad_map_flag


def run_tasks(task):
    index, rxn = task
    return index, drop_unmapped_reactants(rxn)


if __name__ == '__main__':

    # 移除未映射分子
    cleaned_uspto_dir = os.path.join('data', 'Cleaned_data')
    if not os.path.exists(cleaned_uspto_dir):
        os.makedirs(cleaned_uspto_dir, exist_ok=True)
    uspto_remapped_dir = os.path.join('data', 'USPTO_remapped')
    file_list = os.listdir(uspto_remapped_dir)
    all_reactions_df = pd.DataFrame({})
    for file_name in file_list:
        file_path = os.path.join(uspto_remapped_dir, file_name)
        cur_df = pd.read_csv(file_path, sep='\t')
        print(f'读取 {file_name}，共 {len(cur_df)} 条反应')
        all_reactions_df = all_reactions_df.append(cur_df)
    print(f'合并所有反应，总计 {len(all_reactions_df)} 条。')

    mapped_rxn_list = all_reactions_df['mapped_rxn'].tolist()

    tasks = [(i, rxn) for i, rxn in enumerate(mapped_rxn_list)]
    if debug:
        print('调试模式，仅处理前1000条数据……')
        tasks = tasks[:1000]
    run_results = Parallel(n_jobs=6, verbose=1)(
        delayed(run_tasks)(task) for task in tqdm(tasks, desc='处理反应'))

    run_results.sort(key=lambda x: x[0])
    drop_unmapped_rxns, clean_map_rxns = [], []
    bad_map_count = 0

    for idx, (drop_unmapped_rxn, clean_map_rxn, bad_map_flag) in run_results:
        if not bad_map_flag:
            drop_unmapped_rxns.append(drop_unmapped_rxn)
            clean_map_rxns.append(clean_map_rxn)
        else:
            drop_unmapped_rxns.append('')
            clean_map_rxns.append('')
            bad_map_count += 1
            print(f'第{idx}条: {mapped_rxn_list[idx]}，反应物与产物的映射不一致，已丢弃')
    print(f'不一致映射反应总数: {bad_map_count}')
    if not debug:
        all_reactions_df['drop_unmapped_rxn'] = drop_unmapped_rxns
        all_reactions_df['clean_map_rxn'] = clean_map_rxns
        all_reactions_df.to_csv(os.path.join(cleaned_uspto_dir, 'all_clean_remapped_USPTO.csv'), index=False)

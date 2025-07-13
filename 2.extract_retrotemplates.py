import csv
import os

import pandas as pd
# from joblib import Parallel, delayed
import signal
from collections import defaultdict
import multiprocessing
from rdchiral.template_extractor import extract_from_reaction
from tqdm import tqdm

from utils import timeout, get_writer

debug = False





def get_tpl(task):
    clean_map_rxn, patent_id, rxn_smiles = task
    react, reagent, prod = rxn_smiles.split('>')
    reaction = {'_id': patent_id, 'reactants': react, 'products': prod}
    try:
        with timeout(seconds=6):
            template = extract_from_reaction(reaction)
    except Exception as e:
        print(f'提取模板出错: {e}，专利ID: {patent_id}')
        template = None
    return clean_map_rxn, patent_id, rxn_smiles, template


if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    if debug:
        database = pd.read_csv(os.path.join('data', 'Cleaned_data', 'all_clean_remapped_USPTO.csv'), nrows=10000)
    else:
        database = pd.read_csv(os.path.join('data', 'Cleaned_data', 'all_clean_remapped_USPTO.csv'))

    drop_unmapped_rxns = database['drop_unmapped_rxn'].tolist()
    clean_map_rxns = database['clean_map_rxn'].tolist()
    patent_ids = database['source'].tolist()

    unique_rxn_dict = defaultdict(list)
    for patent_id, rxn, clean_map_rxn in zip(patent_ids, drop_unmapped_rxns, clean_map_rxns):
        if not pd.isna(clean_map_rxn):
            if clean_map_rxn not in unique_rxn_dict:
                unique_rxn_dict[clean_map_rxn] = [patent_id, rxn]
    print(f'共计 {len(unique_rxn_dict)} 条唯一反应')

    tasks = []
    for clean_map_rxn, (patent_id, rxn) in list(unique_rxn_dict.items()):
        tasks.append((clean_map_rxn, patent_id, rxn))

    # run_results = Parallel(n_jobs=6, verbose=1)(
    #     delayed(get_tpl)(task) for task in tqdm(tasks))
    fout, writer = get_writer(os.path.join('data', 'Cleaned_data', 'USPTO_remapped_rxn_templates.csv'),
                              ['source', 'droped_unmapped_rxn', 'retro_template', 'clean_map_rxn'])
    # run_results = []
    bad_templates = 0
    good_templates = 0
    with timeout(seconds=10800, error_message='长时间无响应，强制退出。'):
        try:
            for result in tqdm(pool.imap_unordered(get_tpl, tasks), total=len(tasks), desc='提取模板'):
                # run_results.append(result)
                clean_map_rxn, patent_id, rxn, template = result
                if template:
                    if 'reaction_smarts' in template:
                        writer.writerow([patent_id, rxn, template['reaction_smarts'], clean_map_rxn])
                        fout.flush()
                        good_templates += 1
                    else:
                        bad_templates += 1
                else:
                    bad_templates += 1
        except Exception as e:
            print(f'批量提取模板异常: {e}')
    fout.close()
    print(f'模板提取失败 {bad_templates} 条，成功 {good_templates} 条')

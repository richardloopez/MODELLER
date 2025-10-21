#!/usr/bin/env python3.8
# homology_modeling.py

import os
from typing import List, Tuple, Dict, Any

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.selection import Selection
from modeller.parallel import Job, LocalWorker

import config
from config import ALIGN_CODE_TEMPLATE, ALIGN_CODE_SEQUENCE, NUM_MODELS_AUTO, NUM_MODELS_TO_REFINE

def run_automodel(env: Environ, align_file: str, job: Job) -> List[str]:
    """
    Ejecuta AutoModel, genera los modelos base, los renombra y retorna
    los nombres de los Top N modelos seleccionados por DOPEHR.
    """
    initial_models_for_loop_names: List[str] = []
    
    print(f"\n[STEP 4.1] Iniciando AutoModel (Relleno de Gaps) con {NUM_MODELS_AUTO} modelos...")
    
    a = AutoModel(env, 
                  alnfile=align_file, 
                  knowns=ALIGN_CODE_TEMPLATE, 
                  sequence=ALIGN_CODE_SEQUENCE, 
                  assess_methods=(assess.DOPEHR, assess.GA341))
    
    a.use_parallel_job(job)
    a.starting_model = 1
    a.ending_model = NUM_MODELS_AUTO 
    a.library_schedule = autosched.slow
    a.max_var_iterations = 1000
    a.make()
    
    results_auto = a.outputs
    if not results_auto: 
        print("[ERROR] AutoModel falló.")
        return []
    
    sorted_auto_models = sorted(results_auto, key=lambda x: x['DOPE-HR score'])
    
    print(f"\n[STEP 4.1.1] Renombrando los {len(sorted_auto_models)} modelos de AutoModel.")
    
    for model_rank, model_info in enumerate(sorted_auto_models):
        old_name = model_info['name']
        new_name = f'AUTO_{model_rank+1}.pdb'
        
        try:
            os.rename(old_name, new_name)
            if (model_rank + 1) <= NUM_MODELS_TO_REFINE:
                initial_models_for_loop_names.append(new_name)
        except Exception as e:
            print(f"[ERROR] No se pudo renombrar {old_name} a {new_name}. Error: {e}")
            
    
    print(f"\n[STEP 5] AutoModel completado. {len(initial_models_for_loop_names)} modelos (Top {NUM_MODELS_TO_REFINE}) seleccionados para refinamiento de loops.")
    
    return initial_models_for_loop_names

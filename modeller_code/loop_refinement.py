#!/usr/bin/env python3
# loop_refinement.py

import os
from typing import List, Tuple

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.selection import Selection
from modeller.parallel import Job, LocalWorker

import config
from config import ALIGN_CODE_SEQUENCE, CHAIN_ID, NUM_MODELS_LOOP
from custom_models import *


def run_loop_refinement(env: Environ, job: Job, initial_models_names: List[str], loop_ranges: List[Tuple[int, int]]):
    """
    Ejecuta el refinamiento secuencial de loops con DOPEHR para los modelos base.
    """
    
    if not loop_ranges:
        print("\n[STEP 5.1] Saltando refinamiento de loops: No hay loops flexibles definidos.")
        return
    
    valid_loop_ranges = [r for r in loop_ranges if 4 <= (r[1] - r[0] + 1) <= 30]

    if not valid_loop_ranges:
        print("\n[STEP 5.1] Saltando refinamiento de loops: Ningún loop detectado cumple con la longitud requerida (4-30 residuos).")
        return
    
    print(f"\n[STEP 5.2] Iniciando refinamiento dirigido para {len(valid_loop_ranges)} segmentos válidos...")
    
    for model_index, initial_pdb_file in enumerate(initial_models_names):
        
        base_name = initial_pdb_file.replace('.pdb', '')
        
        print(f"\n --- Procesando Modelo Base #{model_index+1}: {initial_pdb_file} ({base_name}) ---")
        
        current_best_pdb_for_thread = initial_pdb_file
        current_base_name_for_refinment = base_name
        
        # Sequentially refine loops
        for j, (start, end) in enumerate(valid_loop_ranges):
            
            print(f"  > Refinando Loop {j+1}/{len(valid_loop_ranges)}: Residuos {start} a {end}")
            
            try:
                ml = DynamicLoopRefiner(env,
                                        inimodel=current_best_pdb_for_thread,
                                        sequence=ALIGN_CODE_SEQUENCE,
                                        loop_start=start,
                                        loop_end=end,
                                        chain_id=CHAIN_ID)

                ml.use_parallel_job(job)
                ml.loop.starting_model = 1
                ml.loop.ending_model = NUM_MODELS_LOOP
                ml.loop.md_level = refine.slow_large
                ml.loop.assess_methods = (assess.DOPEHR, assess.GA341)
                ml.max_var_iterations = 1000

                ml.make()
                
                loop_models_of_this_step = ml.loop.outputs
                
                if loop_models_of_this_step:
                    sorted_loop_outputs_by_loop_dopeHR = sorted(loop_models_of_this_step, key=lambda x: x.get('DOPE-HR score', 9999999.0))
                    
                    for m, model_info in enumerate(sorted_loop_outputs_by_loop_dopeHR):
                        old_name = model_info['name']
                        new_loop_name = f'{current_base_name_for_refinment}_LOOP{j+1}_R{m+1}.pdb'
                        
                        try:
                            os.rename(old_name, new_loop_name)
                        except Exception as e:
                            print(f"    [ERROR] No se pudo renombrar {old_name} a {new_loop_name}. Error: {e}")
                            
                    best_of_this_loop = sorted_loop_outputs_by_loop_dopeHR[0] 
                    current_best_pdb_for_thread = f'{current_base_name_for_refinment}_LOOP{j+1}_R1.pdb'
                    current_base_name_for_refinment = f'{current_base_name_for_refinment}_LOOP{j+1}_R1' 
                
                else:
                    print(f"    -> Advertencia: DOPEHRLoopModel no generó resultados válidos para {start}-{end}. Usando el modelo inicial anterior.")

            except Exception as e:
                print(f"     > ERROR FATAL en DOPEHRLoopModel para {start}-{end} (Modelo Base #{model_index+1}): {e}")
                continue
            
        print(f"\n[STEP 5.2] Refinamiento de Loops completado para el modelo base #{model_index + 1}.")

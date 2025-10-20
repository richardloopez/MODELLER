# controller.py

import sys
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.selection import Selection
from modeller.parallel import Job, LocalWorker

# Importar configuraciones
import config
from config import ALIGNMENT_FILE, ALIGNMENT_CDE_FILE, NUM_PROCESSORS, USE_MANUAL_ALIGNMENT

# Importar módulos de lógica
import utils
import homology_modeling
import loop_refinement

def main_workflow():
    """Ejecuta el pipeline completo de modelado de Modeller."""
    
    # 1. Configuración de Modeller
    env = Environ()
    env.io.atom_files_directory = ['.', '../atom_files'] 
    
    try:
        log.verbose()
    except Exception as e:
        print(f"[WARNING] Fallo al habilitar Modeller logging. La ejecución continuará. Error: {e}")

    # Configurar el paralelismo
    env.jobs = config.NUM_PROCESSORS
    job = Job()
    
    print(f"[PARALLEL] Configurando {NUM_PROCESSORS} workers locales.")
    for _ in range(NUM_PROCESSORS):
        job.append(LocalWorker()) 
    job.start() 

    # 2. Preparación de Alineamiento
    try:
        cde_line, aligned_template_seq, aligned_target_seq = utils.generate_pir_files(
            env, ALIGNMENT_FILE, ALIGNMENT_CDE_FILE, manual_mode=USE_MANUAL_ALIGNMENT
        )
    except Exception as e:
        print(f"\n[ERROR FATAL] Fallo al generar los archivos PIR. Terminando. Error: {e}")
        sys.exit(1)

    # 3. Detección de Loops
    if not aligned_template_seq or not aligned_target_seq:
        print("\n[ERROR] No se pudo obtener el alineamiento. Terminando la ejecución.")
        sys.exit(1)

    loop_ranges_to_refine = utils.find_missing_residues(aligned_template_seq, aligned_target_seq)

    # 4. Filtrado de Loops Flexibles
    if cde_line:
        # La función get_flexible_missing_ranges ya tiene acceso a las constantes SS2_FILE y sequence_full
        loop_ranges_to_refine = utils.get_flexible_missing_ranges(loop_ranges_to_refine)

    # 5. Modelado por Homología (AutoModel)
    initial_models_names = homology_modeling.run_automodel(env, ALIGNMENT_FILE, job)

    # 6. Refinamiento de Loops
    if initial_models_names:
        loop_refinement.run_loop_refinement(env, job, initial_models_names, loop_ranges_to_refine)

    # Asegurar que todos los procesos paralelos han terminado antes de la evaluación final
    print("[PARALLEL] Todos los procesos de Modeller han finalizado.")

    # 7. Evaluación Final y Ranking
    final_ranking, best_final_model = utils.final_evaluation_and_ranking(env)

    if best_final_model:
        print(f"\nEl modelo de más alta calidad (DOPEHR más negativo) fue: {best_final_model['name']} con un Z-score de {best_final_model['DOPEHR Z-score']:.3f}")


if __name__ == '__main__':
    main_workflow()

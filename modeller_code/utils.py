# utils.py

import os
import re
import csv
from typing import List, Tuple, Dict, Any

# Importaciones necesarias de Modeller
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.selection import Selection
from modeller.parallel import Job, LocalWorker

# Importar constantes del archivo de configuración
import config
from config import SS2_FILE, sequence_full, PDB_TEMPLATE_FILE, ALIGN_CODE_TEMPLATE, ALIGN_CODE_SEQUENCE, pdb_aa

# =================================================================
# UTILIDADES DE ALINEAMIENTO Y PIR
# =================================================================

def extract_ss_from_ss2(ss2_file: str, seq_full: str) -> str:
    """Extrae la estructura secundaria predicha de un archivo PSIPRED SS2."""
    # Lógica idéntica a la del código original
    ss_string = ""
    seq_length = len(seq_full)
    try:
        with open(ss2_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.split()
                if len(parts) >= 3:
                    ss_string += parts[2]
                    if len(ss_string) >= seq_length: break
        
        if len(ss_string) < seq_length:
            print(f"[WARNING] Longitud de SS extraída ({len(ss_string)}) es menor que la longitud de la secuencia objetivo ({seq_length}). Rellenando con 'C'.")
            ss_string += 'C' * (seq_length - len(ss_string))

        ss_string_sliced = ss_string[:seq_length]
        print(f"[STEP 2.1] Estructura Secundaria extraída. Longitud utilizada: {len(ss_string_sliced)} para la línea CDE.")
        return ss_string_sliced
        
    except FileNotFoundError:
        print(f"ERROR: Archivo PSIPRED SS2 '{ss2_file}' no encontrado.")
        return ""

def read_sequences_from_ali_temp(ali_file: str) -> Tuple[str, str]:
    """Lee las secuencias alineadas (incluyendo gaps) de un archivo PIR/ALI."""
    # Lógica idéntica a la del código original
    sequences_raw = ["", ""]
    seq_index = -1
    allowed_chars_re = re.compile(r'[A-Z\-\*]')
    
    try:
        with open(ali_file, 'r') as f:
            lines = f.readlines()
            in_sequence_block = False
            current_raw_sequence = ""
            for line in lines:
                line_stripped = line.strip()
                if line_stripped.startswith('>P1;'):
                    if in_sequence_block:
                        cleaned_sequence = "".join(allowed_chars_re.findall(current_raw_sequence.upper()))
                        sequences_raw[seq_index] = cleaned_sequence.split('*')[0]
                        current_raw_sequence = "" 
                    seq_index += 1
                    in_sequence_block = True
                    continue 
                if in_sequence_block and line_stripped:
                    if line_stripped.startswith(('structureX', 'sequence', 'CDE:', '#')): continue
                    current_raw_sequence += line_stripped 
                    if '*' in line_stripped and seq_index == 1:
                        cleaned_sequence = "".join(allowed_chars_re.findall(current_raw_sequence.upper()))
                        sequences_raw[seq_index] = cleaned_sequence.split('*')[0]
                        break 
        
        if len(sequences_raw[0]) == 0 or len(sequences_raw[1]) == 0:
            raise ValueError("No se pudieron extraer dos secuencias alineadas válidas del archivo PIR.")
        
        return sequences_raw[0], sequences_raw[1]
    except Exception as e:
        raise IOError(f"Error al leer el archivo de alineamiento PIR. Revise el formato. Error: {e}")

def generate_pir_files(env: Environ, align_file_modeller: str, align_file_cde: str, manual_mode: bool) -> Tuple[str, str, str]:
    """Genera los archivos PIR finales (con y sin línea CDE) necesarios para Modeller."""
    
    aligned_template_seq = ""
    aligned_target_seq = ""
    cde_line_full = ""
    
    if manual_mode:
        # --- MODO MANUAL: Copiar y leer los archivos PIR proporcionados ---
        print(f"\n[STEP 2] Usando Alineamiento Manual: {config.MANUAL_ALIGNMENT_FILE}")
        # Lógica idéntica a la del código original, usando 'config'
        try:
            os.system(f'cp {config.MANUAL_ALIGNMENT_FILE} {align_file_modeller}')
            aligned_template_seq, aligned_target_seq = read_sequences_from_ali_temp(align_file_modeller)
            
            if not os.path.exists(config.MANUAL_ALIGNMENT_CDE_FILE):
                raise FileNotFoundError(f"Se requiere el archivo PIR con CDE: '{config.MANUAL_ALIGNMENT_CDE_FILE}' en modo manual.")
            os.system(f'cp {config.MANUAL_ALIGNMENT_CDE_FILE} {align_file_cde}')
            
            extract_ss_from_ss2(SS2_FILE, sequence_full) 
            cde_line_full = f"# CDE line copied from {config.MANUAL_ALIGNMENT_CDE_FILE} for reference."
            
        except Exception as e:
            print(f"ERROR FATAL en modo manual. Asegúrese de que ambos archivos existen y tienen el formato PIR correcto. Error: {e}")
            return "", "", ""
            
    else:
        # --- MODO AUTOMÁTICO: Generar usando Modeller.salign() ---
        print("\n[STEP 2] Generando Alineamiento Automático con Modeller.salign()")
        
        # Lógica idéntica a la del código original, usando 'config'
        aln = Alignment(env)
        mdl = Model(env, file=PDB_TEMPLATE_FILE)
        aln.append_model(mdl, align_codes=ALIGN_CODE_TEMPLATE, atom_files=PDB_TEMPLATE_FILE)
        aln.append_sequence(sequence_full) 
        aln[1].code = ALIGN_CODE_SEQUENCE
        aln.salign()
        
        temp_ali_file = "temp_modeller_ali.pir"
        try:
            aln.write(file=temp_ali_file, alignment_format='PIR')
            aligned_template_seq, aligned_target_seq = read_sequences_from_ali_temp(temp_ali_file)
        except Exception as e:
            print(f"ERROR FATAL: Problema al escribir o leer el archivo temporal de Modeller. Error: {e}")
            return "", "", ""
        finally:
            if os.path.exists(temp_ali_file): os.remove(temp_ali_file)

        template_title = f">P1;{ALIGN_CODE_TEMPLATE}"
        template_description = f"structureX:{PDB_TEMPLATE_FILE}:1:{config.CHAIN_ID} :{len(pdb_aa)}:{config.CHAIN_ID}:::-1.00:-1.00" 
        fullseq_title = f">P1;{ALIGN_CODE_SEQUENCE}"
        fullseq_description = f"sequence:{ALIGN_CODE_SEQUENCE}:1: :{len(sequence_full)}: :::-1.00:-1.00"

        ss_string_full = extract_ss_from_ss2(SS2_FILE, sequence_full)
        if not ss_string_full: return "", "", ""

        cde_content = ""
        ss_index = 0
        for char in aligned_target_seq:
            if char == '-': 
                cde_content += '.' 
            else:
                if ss_index < len(ss_string_full):
                    cde_content += ss_string_full[ss_index]
                    ss_index += 1
                else: 
                    cde_content += 'C' 

        cde_line_full = "CDE:" + cde_content
        cde_line_commented = "# " + cde_line_full

        pir_content_modeller = [template_title + "\n", template_description + "\n", f"{aligned_template_seq}*\n", fullseq_title + "\n", fullseq_description + "\n", f"{aligned_target_seq}*\n"]
        pir_content_cde = [template_title + "\n", template_description + "\n", f"{aligned_template_seq}*\n", fullseq_title + "\n", fullseq_description + "\n", cde_line_commented + "\n", f"{aligned_target_seq}*\n"]

        with open(align_file_modeller, 'w') as f: f.writelines(pir_content_modeller)
        with open(align_file_cde, 'w') as f: f.writelines(pir_content_cde)

    print(f"\n[STEP 3] Archivo de Alineamiento PIR (LIMPIO) generado para Modeller: {align_file_modeller}")
    return cde_line_full, aligned_template_seq, aligned_target_seq

# =================================================================
# UTILIDADES DE DETECCIÓN DE LOOPS
# =================================================================

def agrupar_rangos(lista: List[int]) -> List[Tuple[int, int]]:
    """Agrupa una lista de números de residuo individuales en rangos secuenciales."""
    # Lógica idéntica a la del código original
    if not lista: return []
    lista = sorted(lista)
    ranges = []; start = lista[0]; end = lista[0]
    for n in lista[1:]:
        if n == end + 1: end = n
        else: ranges.append((start, end)); start = n; end = n
    ranges.append((start, end)); return ranges

def find_missing_residues(aligned_template_seq: str, aligned_target_seq: str) -> List[Tuple[int, int]]:
    """Identifica los segmentos de inserción (loops) en la secuencia objetivo."""
    # Lógica idéntica a la del código original
    missing_residues: List[int] = []
    target_res_num = 0 

    if len(aligned_template_seq) != len(aligned_target_seq):
        print("\n[ERROR FATAL] Las longitudes de las secuencias alineadas son diferentes. No se puede detectar los loops.")
        return []

    for template_char, target_char in zip(aligned_template_seq, aligned_target_seq):
        
        if target_char != '-':
            target_res_num += 1
        
        if template_char == '-' and target_char.isalpha():
            missing_residues.append(target_res_num)
        
    missing_ranges = agrupar_rangos(missing_residues)
    
    range_strings = [f"[{s}-{e}]" if s != e else f"[{s}]" for s, e in missing_ranges]
    print(f"\n[STEP 1] Rangos de residuos faltantes (loops) detectados: {range_strings}")
    
    return missing_ranges

def get_flexible_missing_ranges(missing_ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Filtra los loops para incluir SOLO aquellos residuos predichos como 'Coil' ('C')."""
    # Lógica idéntica a la del código original
    ss_string_full = extract_ss_from_ss2(SS2_FILE, sequence_full) 
    if not ss_string_full: return missing_ranges
        
    flexible_residues_to_refine = set()
    for start, end in missing_ranges:
        for res_num in range(start, end + 1):
            ss_index = res_num - 1
            if ss_index < len(ss_string_full) and ss_string_full[ss_index] == 'C':
                flexible_residues_to_refine.add(res_num)

    final_flexible_ranges = agrupar_rangos(list(flexible_residues_to_refine))
    range_strings = [f"[{s}-{e}]" if s != e else f"[{s}]" for s, e in final_flexible_ranges]
    print(f"\n[STEP 4.0] Loops flexibles ('C') dentro de regiones faltantes seleccionados para refinamiento: {range_strings}")
    
    return final_flexible_ranges

# =================================================================
# UTILIDADES DE EVALUACIÓN Y REPORTE
# =================================================================

def final_evaluation_and_ranking(env: Environ) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Escanea, evalúa y rankea todos los PDBs generados."""
    
    # Lógica idéntica a la del código original
    print(f"\n{'='*75}\n[STEP 6] INICIANDO EVALUACIÓN FINAL DE TODOS LOS MODELOS PDB\n")
    
    # Buscar archivos PDB generados por el script (AutoModel o DOPEHRLoopModel)
    pdbs_to_calculate_dopeHR = [f for f in os.listdir() if f.endswith(".pdb") and f != PDB_TEMPLATE_FILE and (f.startswith("AUTO_") or "_LOOP" in f) ]
    
    if not pdbs_to_calculate_dopeHR:
        print("[FINAL] No se encontraron archivos PDB generados para evaluar.")
        return [], {}
    
    final_results: List[Dict[str, Any]] = []

    env.io.atom_files_directory = ['.', '../atom_files'] 
    
    for filename in pdbs_to_calculate_dopeHR:
        try:
            mdl = complete_pdb(env, filename)
            atmsel = Selection(mdl.chains[0])
            dopeHR_score = atmsel.assess_dopehr()
            normalized_dopeHR_zscore = mdl.assess_normalized_dopehr() 
            
            final_results.append({
                'name': filename,
                'DOPEHR score': dopeHR_score,
                'DOPEHR Z-score': normalized_dopeHR_zscore
            })
            
            print(f"  -> Evaluado {filename:<40} | DOPEHR: {dopeHR_score:.3f} | Z-score: {normalized_dopeHR_zscore:.3f}")
            
        except Exception as e:
            print(f"  [ERROR] Falló la evaluación de {filename}. Error: {e}")
            final_results.append({
                'name': filename,
                'DOPEHR score': float('inf'), # Valor alto para que vaya al final
                'DOPEHR Z-score': float('inf')
            })
            continue

    final_ranking = sorted(final_results, key=lambda x: x['DOPEHR score'], reverse=False)
    best_final_models = final_ranking[:config.NUM_BEST_FINAL_MODELS]

    print(f"\n{'='*75}\n[CLASIFICACIÓN FINAL] Los {len(final_ranking)} modelos generados (Top {config.NUM_BEST_FINAL_MODELS} mostrados):\n")
    print(f"{'Rank':<5}{'Nombre del Archivo':<45}{'DOPEHR Score':<15}{'Z-score':<10}")
    print(f"{'-'*75}")

    for i, model in enumerate(best_final_models):
        dopeHR_val = model['DOPEHR score']
        zscore_val = model['DOPEHR Z-score']
        
        print(f"{i+1:<5}{model['name']:<45}{dopeHR_val:.3f}{zscore_val:>10.3f}")
            
    print(f"{'='*75}")
    
    write_csv_report(final_ranking)
    
    return final_ranking, final_ranking[0] if final_ranking else {}

def write_csv_report(final_ranking: List[Dict[str, Any]], filename: str = 'Modeling_Report_Final.csv'):
    """Escribe el reporte final de modelos en un archivo CSV."""
    # Lógica idéntica a la del código original
    data_to_write = [{'Nombre': model['name'], 'DOPEHR_Score': model['DOPEHR score'], 'Normalized_DOPEHR_ZScore': model['DOPEHR Z-score']} 
                     for model in final_ranking]
    
    if not data_to_write:
        print("[WARNING] No hay datos para escribir en el archivo CSV.")
        return

    try:
        with open(filename, 'w', newline='') as csvfile:
            fieldnames = ['Nombre', 'DOPEHR_Score', 'Normalized_DOPEHR_ZScore']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data_to_write)
        print(f"\n[FINAL] Reporte CSV generado exitosamente: {filename} con {len(data_to_write)} modelos.")
    except Exception as e:
        print(f"[ERROR] No se pudo escribir el archivo CSV '{filename}'. Error: {e}")

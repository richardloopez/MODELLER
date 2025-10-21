#!/usr/bin/env python3.8
# utils.py

import os
import re
import csv
from typing import List, Tuple, Dict, Any

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.selection import Selection
from modeller.parallel import Job, LocalWorker

import config
from config import SS2_FILE, sequence_full, PDB_TEMPLATE_FILE, ALIGN_CODE_TEMPLATE, ALIGN_CODE_SEQUENCE, pdb_aa, CHAIN_ID

# =================================================================
# UTILIDADES PARA RESIDUOS HETATM / BLK
# =================================================================

def extract_hetatm_residues(pdb_file: str, chain_id: str) -> List[Dict[str, Any]]:
    """
    Extrae información de residuos HETATM del archivo PDB template.
    
    Retorna una lista de diccionarios con información de cada residuo HETATM:
    - resname: nombre del residuo (ej: BLK, HOH, etc.)
    - resnum: número de residuo
    - chain: cadena
    - position_in_sequence: posición relativa en la secuencia (después de qué residuo ATOM)
    """
    hetatm_residues = []
    last_atom_resnum = 0
    seen_hetatm = set()
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    res_chain = line[21:22].strip()
                    if res_chain == chain_id:
                        try:
                            resnum = int(line[22:26].strip())
                            last_atom_resnum = max(last_atom_resnum, resnum)
                        except ValueError:
                            continue
                
                elif line.startswith('HETATM'):
                    res_chain = line[21:22].strip()
                    if res_chain == chain_id:
                        resname = line[17:20].strip()
                        try:
                            resnum = int(line[22:26].strip())
                        except ValueError:
                            continue
                        
                        res_key = (resname, resnum, res_chain)
                        if res_key not in seen_hetatm:
                            seen_hetatm.add(res_key)
                            hetatm_residues.append({
                                'resname': resname,
                                'resnum': resnum,
                                'chain': res_chain,
                                'position_after_atom_resnum': last_atom_resnum
                            })
        
        hetatm_residues.sort(key=lambda x: x['resnum'])
        
        if hetatm_residues:
            print(f"\n[HETATM] Se detectaron {len(hetatm_residues)} residuos HETATM en {pdb_file} (cadena {chain_id}):")
            for het in hetatm_residues:
                print(f"  - {het['resname']} (resnum={het['resnum']}) después del residuo ATOM {het['position_after_atom_resnum']}")
        
        return hetatm_residues
        
    except FileNotFoundError:
        print(f"[WARNING] Archivo PDB '{pdb_file}' no encontrado. No se detectaron residuos HETATM.")
        return []

def insert_blk_in_alignment(aligned_seq: str, hetatm_residues: List[Dict[str, Any]], 
                            template_pdb_length: int) -> str:
    """
    Inserta caracteres '.' (BLK) en la secuencia alineada para representar residuos HETATM.
    
    Los residuos HETATM se insertan al final de la secuencia alineada, ya que típicamente
    están al final del archivo PDB después de todos los residuos ATOM.
    """
    if not hetatm_residues:
        return aligned_seq
    
    # Contar cuántos residuos no-gap tenemos en la secuencia alineada del template
    residue_count = sum(1 for c in aligned_seq if c != '-')
    
    # Agregar un '.' por cada residuo HETATM al final de la secuencia
    blk_chars = '.' * len(hetatm_residues)
    aligned_seq_with_blk = aligned_seq + blk_chars
    
    print(f"\n[HETATM] Insertando {len(hetatm_residues)} caracteres BLK ('.') al final del alineamiento del template")
    
    return aligned_seq_with_blk

# =================================================================
# UTILIDADES DE ALINEAMIENTO Y PIR
# =================================================================

def extract_ss_from_ss2(ss2_file: str, seq_full: str) -> str:
    """Extrae la estructura secundaria predicha de un archivo PSIPRED SS2."""
    ss_string = ""
    seq_length = len(seq_full)
    try:
        with open(ss2_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): 
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    ss_string += parts[2]
                    if len(ss_string) >= seq_length: 
                        break
        
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
    sequences_raw = ["", ""]
    seq_index = -1
    allowed_chars_re = re.compile(r'[A-Z\-\*\.]')  # Agregado '.' para BLK
    
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
                    if line_stripped.startswith(('structureX', 'sequence', 'CDE:', '#')): 
                        continue
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
    """
    Genera los archivos PIR finales (con y sin línea CDE) necesarios para Modeller.
    Incluye soporte para residuos HETATM (BLK).
    """
    
    aligned_template_seq = ""
    aligned_target_seq = ""
    cde_line_full = ""
    
    if manual_mode:
        print(f"\n[STEP 2] Usando Alineamiento Manual: {config.MANUAL_ALIGNMENT_FILE}")
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
        print("\n[STEP 2] Generando Alineamiento Automático con Modeller.salign()")
        
        # 1. Detectar residuos HETATM en el template para información
        hetatm_residues = extract_hetatm_residues(PDB_TEMPLATE_FILE, CHAIN_ID)
        
        # 2. Generar alineamiento básico con Modeller
        # IMPORTANTE: Cuando env.io.hetatm=True, Modeller YA incluye los HETATM como BLK (.)
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
            if os.path.exists(temp_ali_file): 
                os.remove(temp_ali_file)

        # 3. Verificar si Modeller ya incluyó los BLK automáticamente
        blk_count_in_template = aligned_template_seq.count('.')
        
        if blk_count_in_template > 0:
            print(f"\n[HETATM] Modeller detectó automáticamente {blk_count_in_template} residuos BLK en el template")
            print(f"[HETATM] Estos ya están incluidos en el alineamiento como caracteres '.'")
            
            # Modeller ya incluyó los BLK, usamos las secuencias tal cual
            aligned_template_seq_with_blk = aligned_template_seq
            aligned_target_seq_with_blk = aligned_target_seq
            
            # Ajustar el target para tener gaps donde el template tiene BLK
            # Verificar si necesitamos agregar gaps al target
            if len(aligned_template_seq) != len(aligned_target_seq):
                print(f"[WARNING] Longitudes diferentes: template={len(aligned_template_seq)}, target={len(aligned_target_seq)}")
                # Agregar gaps faltantes al target
                len_diff = len(aligned_template_seq) - len(aligned_target_seq)
                if len_diff > 0:
                    aligned_target_seq_with_blk = aligned_target_seq + ('-' * len_diff)
                    print(f"[HETATM] Agregando {len_diff} gaps al target para igualar longitudes")
        else:
            print(f"\n[HETATM] No se detectaron BLK en el alineamiento automático de Modeller")
            print(f"[HETATM] Esto puede ocurrir si el PDB no tiene HETATM o si env.io.hetatm no está configurado")
            aligned_template_seq_with_blk = aligned_template_seq
            aligned_target_seq_with_blk = aligned_target_seq

        # 4. Calcular la longitud real del template (incluyendo HETATM si están presentes)
        # Contar residuos no-gap en la secuencia del template
        template_residue_count = sum(1 for c in aligned_template_seq_with_blk if c not in ['-', ' '])
        template_total_length = template_residue_count
        
        template_title = f">P1;{ALIGN_CODE_TEMPLATE}"
        template_description = f"structureX:{PDB_TEMPLATE_FILE}:1:{CHAIN_ID}:{template_total_length}:{CHAIN_ID}:::-1.00:-1.00"
        fullseq_title = f">P1;{ALIGN_CODE_SEQUENCE}"
        fullseq_description = f"sequence:{ALIGN_CODE_SEQUENCE}:1::{len(sequence_full)}::::-1.00:-1.00"

        # 6. Generar línea CDE con estructura secundaria
        ss_string_full = extract_ss_from_ss2(SS2_FILE, sequence_full)
        if not ss_string_full: 
            return "", "", ""

        cde_content = ""
        ss_index = 0
        for char in aligned_target_seq_with_blk:
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

        # 7. Escribir archivos PIR
        pir_content_modeller = [
            template_title + "\n", 
            template_description + "\n", 
            f"{aligned_template_seq_with_blk}*\n", 
            fullseq_title + "\n", 
            fullseq_description + "\n", 
            f"{aligned_target_seq_with_blk}*\n"
        ]
        pir_content_cde = [
            template_title + "\n", 
            template_description + "\n", 
            f"{aligned_template_seq_with_blk}*\n", 
            fullseq_title + "\n", 
            fullseq_description + "\n", 
            cde_line_commented + "\n", 
            f"{aligned_target_seq_with_blk}*\n"
        ]

        with open(align_file_modeller, 'w') as f: 
            f.writelines(pir_content_modeller)
        with open(align_file_cde, 'w') as f: 
            f.writelines(pir_content_cde)
        
        # Usar las secuencias actualizadas con BLK para el retorno
        aligned_template_seq = aligned_template_seq_with_blk
        aligned_target_seq = aligned_target_seq_with_blk

    print(f"\n[STEP 3] Archivo de Alineamiento PIR (LIMPIO) generado para Modeller: {align_file_modeller}")
    return cde_line_full, aligned_template_seq, aligned_target_seq

# =================================================================
# UTILIDADES DE DETECCIÓN DE LOOPS
# =================================================================

def agrupar_rangos(lista: List[int]) -> List[Tuple[int, int]]:
    """Agrupa una lista de números de residuo individuales en rangos secuenciales."""
    if not lista: 
        return []
    lista = sorted(lista)
    ranges = []
    start = lista[0]
    end = lista[0]
    for n in lista[1:]:
        if n == end + 1: 
            end = n
        else: 
            ranges.append((start, end))
            start = n
            end = n
    ranges.append((start, end))
    return ranges

def find_missing_residues(aligned_template_seq: str, aligned_target_seq: str) -> List[Tuple[int, int]]:
    """
    Identifica los segmentos de inserción (loops) en la secuencia objetivo.
    Ignora los residuos BLK (representados por '.') en el análisis de loops.
    """
    missing_residues: List[int] = []
    target_res_num = 0 

    if len(aligned_template_seq) != len(aligned_target_seq):
        print("\n[ERROR FATAL] Las longitudes de las secuencias alineadas son diferentes. No se puede detectar los loops.")
        return []

    for template_char, target_char in zip(aligned_template_seq, aligned_target_seq):
        
        # Ignorar residuos BLK - no son parte de la proteína
        if template_char == '.':
            continue
        
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
    ss_string_full = extract_ss_from_ss2(SS2_FILE, sequence_full) 
    if not ss_string_full: 
        return missing_ranges
        
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
    
    print(f"\n{'='*75}\n[STEP 6] INICIANDO EVALUACIÓN FINAL DE TODOS LOS MODELOS PDB\n")
    
    pdbs_to_calculate_dopeHR = [
        f for f in os.listdir() 
        if f.endswith(".pdb") 
        and f != PDB_TEMPLATE_FILE 
        and (f.startswith("AUTO_") or "_LOOP" in f)
    ]
    
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
            
            print(f"  -> Evaluado {filename:<40} | DOPEHR: {dopeHR_score:.3f} | Z-score: {normalized_dopeHR_zscore:.3f}")
            
        except Exception as e:
            print(f"  [ERROR] Falló la evaluación de {filename}. Error: {e}")
            final_results.append({
                'name': filename,
                'DOPEHR score': float('inf'),
                'DOPEHR Z-score': float('inf')
            })
            continue

    final_ranking = sorted(final_results, key=lambda x: x['DOPEHR score'], reverse=False)
    best_final_models = final_ranking[:config.NUM_BEST_FINAL_MODELS]

    print(f"\n{'='*75}")
    print(f"RANKING FINAL - Top {config.NUM_BEST_FINAL_MODELS} Modelos por DOPEHR Score")
    print(f"{'='*75}\n")
    print(f"{'Rank':<6} {'Nombre del Modelo':<45} {'DOPEHR':<12} {'Z-score':<12}")
    print(f"{'-'*75}")
    
    for rank, model_data in enumerate(best_final_models, start=1):
        print(f"{rank:<6} {model_data['name']:<45} {model_data['DOPEHR score']:<12.3f} {model_data['DOPEHR Z-score']:<12.3f}")
    
    print(f"\n{'='*75}\n")

    csv_filename = "final_models_ranking.csv"
    try:
        with open(csv_filename, 'w', newline='') as csvfile:
            fieldnames = ['Rank', 'Model Name', 'DOPEHR Score', 'DOPEHR Z-score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for rank, model_data in enumerate(best_final_models, start=1):
                writer.writerow({
                    'Rank': rank,
                    'Model Name': model_data['name'],
                    'DOPEHR Score': f"{model_data['DOPEHR score']:.3f}",
                    'DOPEHR Z-score': f"{model_data['DOPEHR Z-score']:.3f}"
                })
        print(f"[FINAL] Ranking exportado a: {csv_filename}\n")
    except Exception as e:
        print(f"[WARNING] No se pudo escribir el archivo CSV de ranking. Error: {e}")
    
    best_model = best_final_models[0] if best_final_models else {}
    return final_ranking, best_model

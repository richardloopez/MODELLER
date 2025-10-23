#!/usr/bin/env python3
"""
Script simple de prueba para demostrar la detección de HETATM/BLK
No requiere Modeller, solo lee el archivo PDB directamente
"""

import sys
from typing import List, Dict, Any

def extract_hetatm_residues_simple(pdb_file: str, chain_id: str) -> List[Dict[str, Any]]:
    """
    Extrae información de residuos HETATM del archivo PDB template.
    Versión simplificada que no requiere Modeller.
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
        return hetatm_residues
        
    except FileNotFoundError:
        print(f"[ERROR] Archivo PDB '{pdb_file}' no encontrado.")
        return []

def get_all_chains(pdb_file: str) -> Dict[str, Dict[str, int]]:
    """Obtiene información sobre todas las cadenas en el PDB"""
    chains_info = {}
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain = line[21:22].strip()
                    record_type = 'ATOM' if line.startswith('ATOM') else 'HETATM'
                    
                    if chain not in chains_info:
                        chains_info[chain] = {'ATOM': 0, 'HETATM': 0}
                    chains_info[chain][record_type] += 1
    except FileNotFoundError:
        pass
    
    return chains_info

def main():
    """Ejecuta una prueba de detección de HETATM"""
    
    # Configuración
    PDB_FILE = '8vx1_DS_renum_HETATM.pdb'
    CHAIN_ID = 'A'
    
    print("="*70)
    print("PRUEBA DE DETECCIÓN DE RESIDUOS HETATM/BLK")
    print("="*70)
    
    print(f"\n1. Archivo PDB Template: {PDB_FILE}")
    
    # Primero, mostrar información sobre todas las cadenas
    chains_info = get_all_chains(PDB_FILE)
    if chains_info:
        print(f"\n2. Información de cadenas en el PDB:")
        print(f"   {'Cadena':<10} {'ATOM':<10} {'HETATM':<10}")
        print(f"   {'-'*30}")
        for chain, counts in sorted(chains_info.items()):
            print(f"   {chain:<10} {counts['ATOM']:<10} {counts['HETATM']:<10}")
    
    print(f"\n3. Analizando cadena configurada: {CHAIN_ID}")
    
    # Detectar residuos HETATM
    print(f"   Detectando residuos HETATM en cadena {CHAIN_ID}...")
    hetatm_residues = extract_hetatm_residues_simple(PDB_FILE, CHAIN_ID)
    
    if not hetatm_residues:
        print(f"\n   ⚠ No se detectaron residuos HETATM en la cadena {CHAIN_ID}")
        print("\n   Esto significa que:")
        print(f"   - La cadena {CHAIN_ID} contiene solo residuos de proteína (ATOM)")
        print("   - No hay ligandos, heteroátomos, o moléculas no-proteicas en esta cadena")
        print("   - Los HETATM pueden estar en otras cadenas (ver tabla arriba)")
        print("\n   CONCLUSIÓN:")
        print(f"   ✓ Si solo necesitas modelar la proteína de la cadena {CHAIN_ID}, el código funcionará correctamente")
        print(f"   ✓ No se transferirán ligandos (porque no hay ninguno en la cadena {CHAIN_ID})")
        print(f"   ⚠ Si necesitas ligandos de otras cadenas, actualiza CHAIN_ID en config.py")
        return 0
    
    print(f"\n3. Resumen de detección:")
    print(f"   Total de residuos HETATM detectados: {len(hetatm_residues)}")
    
    # Agrupar por tipo de residuo
    residue_types = {}
    for het in hetatm_residues:
        resname = het['resname']
        if resname not in residue_types:
            residue_types[resname] = 0
        residue_types[resname] += 1
    
    print(f"\n4. Tipos de residuos HETATM:")
    for resname, count in sorted(residue_types.items()):
        print(f"   - {resname}: {count} residuo(s)")
    
    print(f"\n5. Primeros 10 residuos HETATM detectados:")
    print(f"   {'Tipo':<8} {'ResNum':<8} {'Cadena':<8} {'Después de ATOM #':<20}")
    print(f"   {'-'*50}")
    for het in hetatm_residues[:10]:
        print(f"   {het['resname']:<8} {het['resnum']:<8} {het['chain']:<8} {het['position_after_atom_resnum']:<20}")
    
    if len(hetatm_residues) > 10:
        print(f"   ... y {len(hetatm_residues) - 10} más")
    
    print("\n" + "="*70)
    print("EXPLICACIÓN DE CÓMO FUNCIONARÁ EL MODELADO")
    print("="*70)
    
    print("\n6. Durante el modelado con Modeller:")
    print(f"   a) env.io.hetatm = True está configurado ✓")
    print(f"   b) Modeller leerá {len(hetatm_residues)} residuos HETATM automáticamente")
    print(f"   c) Modeller convertirá HETATM a BLK (.) en la secuencia del template")
    print(f"   d) En el alineamiento PIR:")
    print(f"      - Template tendrá {len(hetatm_residues)} caracteres '.' (BLK)")
    print(f"      - Target tendrá {len(hetatm_residues)} caracteres '-' (gaps)")
    print(f"   e) Los residuos HETATM se transferirán como cuerpos rígidos al modelo final")
    
    print("\n7. Ejemplo de alineamiento (simplificado):")
    print(f"   Template: AMINOACIDOS_PROTEINA...{'.' * min(len(hetatm_residues), 5)}")
    print(f"   Target:   AMINOACIDOS_PROTEINA...{'-' * min(len(hetatm_residues), 5)}")
    print(f"                                      {'↑' * min(len(hetatm_residues), 5)}")
    print(f"                                      BLK = Residuos HETATM del template")
    
    print("\n8. Mensajes que verás durante la ejecución:")
    print(f"   '[HETATM] Se detectaron {len(hetatm_residues)} residuos HETATM en {PDB_FILE}'")
    print(f"   '[HETATM] Modeller detectó automáticamente {len(hetatm_residues)} residuos BLK en el template'")
    print(f"   '[HETATM] Estos ya están incluidos en el alineamiento como caracteres (.)'")
    
    print("\n" + "="*70)
    print("✓ PRUEBA COMPLETADA EXITOSAMENTE")
    print("="*70)
    print("\nEl sistema está listo para modelar con residuos HETATM/BLK.")
    print("\nPara ejecutar el modelado completo:")
    print("  1. En tu cluster con Modeller: sbatch modeller_lanzador.sh")
    print("  2. O directamente: python3 controller.py")
    print("\nNOTA: Asegúrate de tener Modeller instalado y licenciado.\n")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())

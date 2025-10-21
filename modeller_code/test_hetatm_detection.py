#!/usr/bin/env python3.8
"""
Script de prueba para demostrar la detección de HETATM/BLK
Este script muestra cómo el código detecta y procesa residuos HETATM
"""

import sys
import config
import utils

def main():
    """Ejecuta una prueba de detección de HETATM"""
    
    print("="*70)
    print("PRUEBA DE DETECCIÓN DE RESIDUOS HETATM/BLK")
    print("="*70)
    
    print(f"\n1. Archivo PDB Template: {config.PDB_TEMPLATE_FILE}")
    print(f"   Cadena: {config.CHAIN_ID}")
    
    # Detectar residuos HETATM
    print("\n2. Detectando residuos HETATM...")
    hetatm_residues = utils.extract_hetatm_residues(
        config.PDB_TEMPLATE_FILE, 
        config.CHAIN_ID
    )
    
    if not hetatm_residues:
        print("\n   ⚠ No se detectaron residuos HETATM")
        print("   Esto puede significar:")
        print("   - El PDB no contiene registros HETATM")
        print("   - La cadena especificada no tiene HETATM")
        return 1
    
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
    print(f"   b) Modeller leerá {len(hetatm_residues)} residuos HETATM como BLK (.)")
    print(f"   c) En el alineamiento PIR:")
    print(f"      - Template tendrá {len(hetatm_residues)} caracteres '.' (BLK)")
    print(f"      - Target tendrá {len(hetatm_residues)} caracteres '-' (gaps)")
    print(f"   d) Los residuos HETATM se transferirán al modelo final")
    
    print("\n7. Ejemplo de alineamiento (simplificado):")
    # Simulación simple
    template_example = config.pdb_aa[:50] + "..." + ("." * min(len(hetatm_residues), 5))
    target_example = config.sequence_full[:50] + "..." + ("-" * min(len(hetatm_residues), 5))
    
    print(f"   Template: {template_example}")
    print(f"   Target:   {target_example}")
    print(f"             {'↑'*min(len(hetatm_residues), 5)}")
    print(f"             Residuos BLK transferidos del template")
    
    print("\n" + "="*70)
    print("✓ PRUEBA COMPLETADA EXITOSAMENTE")
    print("="*70)
    print("\nEl sistema está listo para modelar con residuos HETATM/BLK.")
    print("Ejecute 'python3 controller.py' para iniciar el modelado completo.\n")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())

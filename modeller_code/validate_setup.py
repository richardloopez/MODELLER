#!/usr/bin/env python3
"""
Script de validación para el pipeline de Modeller con soporte HETATM/BLK
Verifica que todos los componentes están instalados y configurados correctamente.
"""

import sys
import os

def validate_files():
    """Verifica que todos los archivos necesarios existen"""
    print("="*70)
    print("VALIDACIÓN DE ARCHIVOS")
    print("="*70)
    
    required_files = [
        '8vx1_DS_renum_HETATM.pdb',
        'P1_DHX_secondary_strucuture.ss2',
        'controller.py',
        'config.py',
        'utils.py',
        'homology_modeling.py',
        'loop_refinement.py',
        'custom_models.py'
    ]
    
    all_exist = True
    for f in required_files:
        exists = os.path.exists(f)
        status = "✓ OK" if exists else "✗ FALTA"
        print(f"  {status:8} {f}")
        if not exists:
            all_exist = False
    
    print()
    return all_exist

def validate_imports():
    """Verifica que se pueden importar todos los módulos"""
    print("="*70)
    print("VALIDACIÓN DE IMPORTACIONES")
    print("="*70)
    
    try:
        print("  Importando modeller...", end=" ")
        from modeller import Environ, Alignment, Model
        from modeller.automodel import AutoModel, assess, autosched
        print("✓ OK")
    except ImportError as e:
        print(f"✗ ERROR: {e}")
        return False
    
    try:
        print("  Importando módulos del proyecto...", end=" ")
        import config
        import utils
        import homology_modeling
        import loop_refinement
        import custom_models
        print("✓ OK")
    except ImportError as e:
        print(f"✗ ERROR: {e}")
        return False
    
    print()
    return True

def validate_hetatm_detection():
    """Verifica que se pueden detectar residuos HETATM"""
    print("="*70)
    print("VALIDACIÓN DE DETECCIÓN HETATM")
    print("="*70)
    
    try:
        import utils
        import config
        
        print(f"\n  Analizando archivo: {config.PDB_TEMPLATE_FILE}")
        print(f"  Cadena: {config.CHAIN_ID}\n")
        
        hetatm_residues = utils.extract_hetatm_residues(
            config.PDB_TEMPLATE_FILE, 
            config.CHAIN_ID
        )
        
        if hetatm_residues:
            print(f"\n  ✓ Se detectaron {len(hetatm_residues)} residuos HETATM correctamente")
            print("\n  Detalles:")
            for het in hetatm_residues[:10]:  # Mostrar solo los primeros 10
                print(f"    - {het['resname']:5} resnum={het['resnum']:4} (después de ATOM resnum {het['position_after_atom_resnum']})")
            if len(hetatm_residues) > 10:
                print(f"    ... y {len(hetatm_residues) - 10} más")
        else:
            print("  ⚠ No se detectaron residuos HETATM")
        
        print()
        return True
        
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        return False

def validate_environment():
    """Verifica el entorno de Modeller"""
    print("="*70)
    print("VALIDACIÓN DE ENTORNO MODELLER")
    print("="*70)
    
    try:
        from modeller import Environ
        env = Environ()
        env.io.hetatm = True
        print("  ✓ env.io.hetatm = True configurado correctamente")
        
        env.io.atom_files_directory = ['.', '../atom_files']
        print("  ✓ Directorios de archivos configurados")
        
        print()
        return True
        
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        return False

def show_alignment_preview():
    """Muestra una vista previa de cómo se vería el alineamiento con BLK"""
    print("="*70)
    print("VISTA PREVIA DE ALINEAMIENTO CON BLK")
    print("="*70)
    
    try:
        import utils
        import config
        
        # Detectar HETATM
        hetatm_residues = utils.extract_hetatm_residues(
            config.PDB_TEMPLATE_FILE, 
            config.CHAIN_ID
        )
        
        # Simular secuencias alineadas (primeros 60 caracteres como ejemplo)
        template_seq_example = config.pdb_aa[:60].replace('', '')  # Primeros 60 residuos
        target_seq_example = config.sequence_full[:60].replace('', '')
        
        # Agregar BLK al final
        template_with_blk = template_seq_example + ('.' * len(hetatm_residues))
        target_with_gaps = target_seq_example + ('-' * len(hetatm_residues))
        
        print("\n  Ejemplo de alineamiento (primeros 60 residuos + BLK):\n")
        print(f"  Template: {template_with_blk[:70]}")
        print(f"  Target:   {target_with_gaps[:70]}")
        
        if len(hetatm_residues) > 0:
            print(f"\n  Los últimos {len(hetatm_residues)} caracteres '.' representan residuos HETATM/BLK")
            print(f"  Los últimos {len(hetatm_residues)} caracteres '-' en target son gaps para alinear")
        
        print()
        return True
        
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        return False

def main():
    """Ejecuta todas las validaciones"""
    print("\n" + "="*70)
    print("VALIDACIÓN DEL PIPELINE DE MODELLER CON SOPORTE HETATM/BLK")
    print("="*70 + "\n")
    
    results = []
    
    results.append(("Archivos", validate_files()))
    results.append(("Importaciones", validate_imports()))
    results.append(("Entorno Modeller", validate_environment()))
    results.append(("Detección HETATM", validate_hetatm_detection()))
    results.append(("Vista previa alineamiento", show_alignment_preview()))
    
    # Resumen final
    print("="*70)
    print("RESUMEN DE VALIDACIÓN")
    print("="*70)
    
    all_passed = True
    for name, passed in results:
        status = "✓ PASADO" if passed else "✗ FALLADO"
        print(f"  {status:12} {name}")
        if not passed:
            all_passed = False
    
    print("\n" + "="*70)
    if all_passed:
        print("✓ TODAS LAS VALIDACIONES PASARON EXITOSAMENTE")
        print("\nEl pipeline está listo para ejecutarse.")
        print("\nPara ejecutar el modelado completo:")
        print("  - En un clúster con SLURM: sbatch modeller_lanzador.sh")
        print("  - Localmente: python3 controller.py")
    else:
        print("✗ ALGUNAS VALIDACIONES FALLARON")
        print("\nPor favor revise los errores arriba.")
    print("="*70 + "\n")
    
    return 0 if all_passed else 1

if __name__ == '__main__':
    sys.exit(main())

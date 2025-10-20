#!/usr/bin/env python3.8
# custom_models.py

from modeller import *
from modeller.automodel import DOPEHRLoopModel # Clase base para refinamiento de loops
from modeller.selection import Selection

# Definición de la clase en el ámbito GLOBAL
class DynamicLoopRefiner(DOPEHRLoopModel):
    """
    Clase personalizada de DOPEHRLoopModel, definida globalmente
    para permitir el procesamiento paralelo (pickling).
    """

    def __init__(self, env, inimodel, sequence, loop_start, loop_end, chain_id, **kwargs):
        # La clase hereda de DOPEHRLoopModel
        super().__init__(env,
                         inimodel=inimodel,
                         sequence=sequence,
                         **kwargs)
        # Almacenamos los límites del loop y la cadena en la instancia (self.)
        self.loop_start = loop_start
        self.loop_end = loop_end
        self.chain_id = chain_id

    def select_loop_atoms(self):
        """Define los residuos que serán refinados usando los atributos de la instancia."""

        range_start = f'{self.loop_start}:{self.chain_id}'
        range_end = f'{self.loop_end}:{self.chain_id}'
        
        # residue_range() selecciona un rango específico de residuos
        return Selection(self.residue_range(range_start, range_end))

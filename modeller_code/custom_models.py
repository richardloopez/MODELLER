#!/usr/bin/env python3.8
# custom_models.py

from modeller import *
from modeller.automodel import DOPEHRLoopModel
from modeller.selection import Selection

class DynamicLoopRefiner(DOPEHRLoopModel):
    """
    Clase personalizada de DOPEHRLoopModel, definida globalmente
    para permitir el procesamiento paralelo (pickling).
    """

    def __init__(self, env, inimodel, sequence, loop_start, loop_end, chain_id, **kwargs):
        super().__init__(env,
                         inimodel=inimodel,
                         sequence=sequence,
                         **kwargs)
        self.loop_start = loop_start
        self.loop_end = loop_end
        self.chain_id = chain_id

    def select_loop_atoms(self):
        """Define los residuos que serán refinados usando los atributos de la instancia."""
        range_start = f'{self.loop_start}:{self.chain_id}'
        range_end = f'{self.loop_end}:{self.chain_id}'
        return Selection(self.residue_range(range_start, range_end))

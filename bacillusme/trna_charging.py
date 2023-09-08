from __future__ import print_function, absolute_import, division

from six import iteritems

import cobrame
from bacillusme.corrections import correct_trna_modifications

amino_acid_trna_synthetase = {
  "cys__L_c": "BSU00940-MONOMER",
  "leu__L_c": "BSU30320-MONOMER",
  "lys__L_c": "BSU00820-MONOMER",
  "asp__L_c": "BSU27550-MONOMER",
  "phe__L_c": "CPLX8J2-11",
  "his__L_c": "BSU27560-MONOMER",
  "asn__L_c": "BSU22360-MONOMER",
  "pro__L_c": "BSU16570-MONOMER",
  "ala__L_c": "BSU27410-MONOMER",
  "ile__L_c": "BSU15430-MONOMER",
  "ser__L_c": "BSU00130-MONOMER",
  "arg__L_c": "BSU37330-MONOMER",
  "met__L_c": "BSU00380-MONOMER",
  "tyr__L_c": "BSU29670-MONOMER",
  "glu__L_c": "CPLX8J2-4",
  "thr__L_c": "BSU28950-MONOMER",
  "val__L_c": "BSU28090-MONOMER",
  "gly_c": "CPLX8J2-12",
  "trp__L_c": "BSU11420-MONOMER",
  "gln__L_c": "CPLX8J2-4"
}

trna_modification = {}

modification_info = {}


def add_trna_modification_procedures(model):

    modifications = trna_modification.copy()
    modifications = correct_trna_modifications(modifications)

    for mod, components in iteritems(modifications):
        trna_mod = cobrame.SubreactionData(mod, model)
        trna_mod.enzyme = components['machines']
        trna_mod.stoichiometry = components['metabolites']
        trna_mod.keff = 65.  # iOL uses 65 for all tRNA mods
        if 'carriers' in components.keys():
            for carrier, stoich in components['carriers'].items():
                if stoich < 0:
                    trna_mod.enzyme += [carrier]
                trna_mod.stoichiometry[carrier] = stoich

        # Add element contribution from modification to tRNA
        trna_mod._element_contribution = \
            modification_info[mod.split('_')[0]]['elements']

    return modifications

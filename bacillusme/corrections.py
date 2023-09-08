from __future__ import print_function, absolute_import, division

from six import iteritems
from warnings import warn
import json

import pandas as pd

import cobrame

# MetabolicReactions present in iOL1650 that have been removed for iJL1678b
removed_reactions = []

# PFL enzymes are activated by a glycyl radical group
pfl_isozymes = []

making_iOL_model = False


def update_metabolite_formulas(m_model):
    # Formulas of metabolites not included in metabolites.txt
    formulas = [('4fe4s_c', 'Fe4S4'), ('2fe2s_c', 'Fe2S2')]

    for met, formula in formulas:
        try:
            met_obj = m_model.metabolites.get_by_id(met)
        except KeyError:
            warn('Creating new metabolite (%s)' % met)
            met_obj = cobrame.Metabolite(met)
            m_model.add_metabolites([met_obj])
        met_obj.formula = formula


def correct_complex_modifications(model):

    return


def correct_reaction_matrix(reaction_matrix_dict):

    return reaction_matrix_dict


def correct_reaction_info_frame(df):
    return df


def correct_reaction_stoichiometries(model, file_name):
    df = pd.read_excel(file_name, index_col=0)
    for d in df.index:
        stoich_data = model.process_data.get_by_id(d)
        stoich_dict = json.loads(df.loc[d, 'Stoich Change'].replace("'", "\""))
        for met, stoich in iteritems(stoich_dict):
            stoich_data.stoichiometry[met] = stoich

        stoich_data.update_parent_reactions()


def correct_tu_dataframe(df):
    return df


def correct_enzyme_reaction_association_frame(df):
    return df


def correct_trna_modifications(mod):
    return mod


def correct_rrna_modifications(mod):
    if making_iOL_model:
        return mod
    return mod


def correct_complex_stoichiometry(stoichiometry):
    # Corrections to subunit stoichiometry.

    return stoichiometry


def correct_complex_modification_dict(stoichiometry):
    # specific patches.

    return stoichiometry

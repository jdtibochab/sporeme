from __future__ import print_function, absolute_import, division


import cobrame
from bacillusme import corrections

# these guys can transfer assembled iron sulfur clusters to the various enzymes
fes_transfer = {'BSU32680-MONOMER':'BSU32680-MONOMER'}

# Add known specific chaperone
fes_chaperones = {}  # FE-S modification

# complexes that can transfer an iron sulfur cluster to an enzyme target
generic_fes_transfer_complexes = {
    'generic_2fe2s_transfer_complex': ['BSU32680-MONOMER_mod_1:2fe2s'],
    'generic_4fe4s_transfer_complex': ['BSU32680-MONOMER_mod_1:4fe4s']
} #sufU Albrecht et al, 2010

lipoate_modifications = {}

bmocogdp_chaperones = {}


def add_iron_sulfur_modifications(me_model):

    for name, complexes in generic_fes_transfer_complexes.items():
        generic_fes_transfer = cobrame.GenericData(name, me_model, complexes)
        generic_fes_transfer.create_reactions()

    for fes in ['2fe2s_c', '4fe4s_c']:
        me_model.add_metabolites([cobrame.Metabolite(fes)])
        for name in fes_transfer.values():
            rxn = cobrame.MEReaction('_'.join([name, fes, 'unloading']))
            me_model.add_reactions([rxn])
            rxn.add_metabolites({name + '_mod_1:' + fes.replace('_c', ''): -1,
                                 fes: 1,
                                 name: 1})

    # add fes transfer enzymes to proper modification data
    mod_2fe2s = me_model.process_data.mod_2fe2s_c
    mod_2fe2s.enzyme = 'generic_2fe2s_transfer_complex'
    mod_2fe2s.stoichiometry = {'2fe2s_c': -1.}

    mod_4fe4s = me_model.process_data.mod_4fe4s_c
    mod_4fe4s.enzyme = 'generic_4fe4s_transfer_complex'
    mod_4fe4s.stoichiometry = {'4fe4s_c': -1.}

    for chaperone in set(fes_chaperones.values()):
        new_mod = cobrame.SubreactionData('mod_2fe2s_c_' + chaperone, me_model)
        new_mod.enzyme = [chaperone, 'generic_2fe2s_transfer_complex']
        new_mod.stoichiometry = {'2fe2s_c': -1.}

    for cplx_data in me_model.process_data.get_by_id(
            'mod_2fe2s_c').get_complex_data():
        cplx_id = cplx_data.id.split('_mod')[0]
        if cplx_id in fes_chaperones:
            cplx_data.subreactions['mod_2fe2s_c_' + fes_chaperones[
                cplx_id]] = \
                cplx_data.subreactions.pop('mod_2fe2s_c')


def add_lipoate_modifications(me_model):
    return


def add_bmocogdp_modifications(me_model):
    return


def add_modification_procedures(me_model):
    # add SubreactionData for iron sulfur clusters
    add_iron_sulfur_modifications(me_model)

    # lipoate modifications can be accomplished using two different mechanisms
    add_lipoate_modifications(me_model)

    # bmocogdp modifications have multiple selective chaperones that transfer
    # the metabolite to the target complexes
    add_bmocogdp_modifications(me_model)

    corrections.correct_complex_modifications(me_model)

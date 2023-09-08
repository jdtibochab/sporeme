from __future__ import print_function, absolute_import, division

from six import iteritems

from cobrame import ComplexData, SubreactionData
from bacillusme.corrections import correct_rrna_modifications


# Dictionary of {formation_step:[{metabolite:stoichiometry}]}
# Positive for reactants negative for products (complex formation convention)

# Leaving out Tig_mono trigger factor for now. It's a chaperone not part of
# formation

###### From Teddy's notes in original ME code #######
# What I'm trying to accomplish here: (old, but still useful)
# 1) 1 30Sp + 1 generic_16S -> 1 rib_30
# 2) 1 50Sp + 1 generic_23S + 1 generic_5S -> 1 rib_50
# 3) 1 rib_30 + 1 rib_50 -> 1 rib_70
# 4) 1 rib_70 + 1 b0884_assumedMonomer + 1 b1718_uniprotComplex  -->
#       1 rib_30_IF1_IF3 + 1 rib_50
# 5) 1 b3168_assumedMonomer_gtp + 1 rib_30_IF1_IF3 --> 1 rib_30_ini


# Mod:
#   Era_dim (assembly factor) + 2 gtp +2 h20->
#   30S_assembly_factor_gtp_hydrolying_assembly_phase_1_gtp ->
#           2 gtp + 2 h + 2 pi


# TODO Check how 2 gtp is added
# Ribosome stoichiometry from SubtiWiki AND BsubCyc
ribosome_stoich = {'30_S_assembly': {'stoich': {'BSU22880-MONOMER': 1, # S1
                                                'BSU16490-MONOMER': 1, # S2
                                                'BSU01220-MONOMER': 1, # S3
                                                'BSU29660-MONOMER': 1, # S4
                                                'BSU01330-MONOMER': 1, # S5
                                                'BSU40910-MONOMER': 1, # S6
                                                'BSU01110-MONOMER': 1, # S7
                                                'BSU01300-MONOMER': 1, # S8
                                                'BSU01500-MONOMER': 1, # S9
                                                'BSU01150-MONOMER': 1, # S10
                                                'BSU01420-MONOMER': 1, # S11
                                                'BSU01100-MONOMER': 1, # S12
                                                'BSU01410-MONOMER': 1, # S13
                                                'BSU01290-MONOMER': 1, # S14
                                                'BSU16680-MONOMER': 1, # S15
                                                'BSU15990-MONOMER': 1, # S16
                                                'BSU01250-MONOMER': 1, # S17
                                                'BSU40890-MONOMER': 1, # S18
                                                'BSU01200-MONOMER': 1, # S19
                                                'BSU25550-MONOMER': 1, # S20
                                                'BSU25410-MONOMER': 1, # S21                      
                                                'generic_16s_rRNAs':1}},
                    '50_S_assembly': {'stoich': {'generic_23s_rRNAs': 1,
                                                  'generic_5s_rRNAs': 1,
                                                  'BSU01030-MONOMER': 1, # L1
                                                  'BSU01190-MONOMER': 1, # L2
                                                  'BSU01160-MONOMER': 1, # L3
                                                  'BSU01170-MONOMER': 1, # L4
                                                  'BSU01280-MONOMER': 1, # L5
                                                  'BSU01310-MONOMER': 1, # L6
                                                  'BSU40500-MONOMER': 1, # L9
                                                  'BSU01040-MONOMER': 1, # L10
                                                  'BSU01020-MONOMER': 1, # L11
                                                  'BSU01050-MONOMER': 1, # L12
                                                  'BSU01490-MONOMER': 1, # L13
                                                  'BSU01260-MONOMER': 1, # L14
                                                  'BSU01350-MONOMER': 1, # L15
                                                  'BSU01230-MONOMER': 1, # L16
                                                  'BSU01440-MONOMER': 1, # L17
                                                  'BSU01320-MONOMER': 1, # L18
                                                  'BSU16040-MONOMER': 1, # L19
                                                  'BSU28850-MONOMER': 1, # L20
                                                  'BSU27960-MONOMER': 1, # L21
                                                  'BSU01210-MONOMER': 1, # L22
                                                  'BSU01180-MONOMER': 1, # L23
                                                  'BSU01270-MONOMER': 1, # L24
                                                  'BSU27940-MONOMER': 1, # L27
                                                  'BSU15820-MONOMER': 1, # L28
                                                  'BSU01240-MONOMER': 1, # L29
                                                  'BSU01340-MONOMER': 1, # L30
                                                  'BSU37070-MONOMER': 1, # L31
                                                  'BSU15080-MONOMER': 1, # L32
                                                  'BSU24900-MONOMER': 1, # L33a
                                                  'BSU00990-MONOMER': 1, # L33b
                                                  'BSU41060-MONOMER': 1, # L34
                                                  'BSU28860-MONOMER': 1, # L35
                                                  'BSU01400-MONOMER': 1, # L36
                                                  'BSU28230-MONOMER': 1}, # Trigger factor
                                                  'mods': None,
                                                  'enzymes': None},
                   # TODO Make sure this isn't double counted
                   'assemble_ribosome_subunits': {'stoich': {'gtp_c': 1}
                                                  }}

ribosome_subreactions = {'gtp_bound_30S_assembly_factor_phase1':
                         {'enzyme': 'BSU25290-MONOMER',
                          'stoich': {'gtp_c': -2,
                                     'h2o_c': -2,
                                     'h_c': 2,
                                     'pi_c': 2,
                                     'gdp_c': 2},
                          'num_mods': 1},

                         'RbfA_mono_assembly_factor_phase1':
                         {'enzyme': 'BSU16650-MONOMER',
                          'stoich': {},
                          'num_mods': 1},

                         'RimM_mono_assembly_factor_phase1':
                         {'enzyme': 'BSU16020-MONOMER',
                          'stoich': {},
                          'num_mods': 1}
                         }


def add_ribosome(me_model, verbose=True):
    ribosome_complex = ComplexData("ribosome", me_model)
    ribosome_components = ribosome_complex.stoichiometry

    rrna_mods = correct_rrna_modifications(rrna_modifications)
    for mod, components in iteritems(rrna_mods):
        rrna_mod = SubreactionData(mod, me_model)
        rrna_mod.enzyme = components['machine']
        rrna_mod.stoichiometry = components['metabolites']
        rrna_mod.keff = 65.  # iOL uses 65. for all RNA mods

        # Add element contribution from modification to rRNA
        rrna_mod._element_contribution = \
            modification_info[mod.split('_')[0]]['elements']

        if 'carriers' in components.keys():
            for carrier, stoich in iteritems(components['carriers']):
                if stoich < 0:
                    rrna_mod.enzyme += [carrier]
                rrna_mod.stoichiometry[carrier] = stoich
        ribosome_complex.subreactions[rrna_mod.id] = 1

    subreaction_dict = ribosome_subreactions
    for subreaction_id in subreaction_dict:
        # get subreaction info
        subreaction_stoich = subreaction_dict[subreaction_id]['stoich']
        subreaction_enzyme = subreaction_dict[subreaction_id]['enzyme']
        num_subreactions = subreaction_dict[subreaction_id]['num_mods']

        # add subreaction to model
        subreaction = SubreactionData(subreaction_id, me_model)
        subreaction.stoichiometry = subreaction_stoich
        subreaction.enzyme = subreaction_enzyme

        # account for subreactions in complex data
        ribosome_complex.subreactions[subreaction.id] = num_subreactions

    # Ribosomes in iOL1650 contain 171 mg2 ions
    ribosome_complex.subreactions['mod_mg2_c'] = 171.
    ribosome_assembly = ribosome_stoich
    for process in ribosome_assembly:
        for protein, amount in iteritems(ribosome_assembly[process]['stoich']):
            ribosome_components[protein] += amount
    ribosome_complex.create_complex_formation(verbose=verbose)

rrna_modifications = {
                      # ---------16S Modifications---------------
                      'm4Cm_at_1402': {'machine': 'generic_16Sm4Cm1402',
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},

                      # ---------23S Modifications---------------
                      'm2A_at_2503': {'machine': 'BSU15750-MONOMER',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5U_at_747': {'machine': 'BSU06730-MONOMER',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}}}

modification_info = {'Y': {'elements': {}, 'charge': 0},
                     'm3Y': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'Um': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'm7G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm6A': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm5U': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm5C': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm2G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm2A': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm1G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'Gm': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'D': {'elements': {'H': 2}, 'charge': 0},
                     'Cm': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm62A': {'elements': {'C': 2, 'H': 4}, 'charge': 0},
                     'm4Cm': {'elements': {'C': 2, 'H': 4}, 'charge': 0},
                     'm3U': {'elements': {'C': 1, 'H': 2}, 'charge': 0}
                     }

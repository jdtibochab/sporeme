from __future__ import division, absolute_import, print_function

from cobrame import SubreactionData, Complex
from cobrame.util import dogma

# 1 machine + 1 atp + 1 aa + 1 h2o --> 1 machine-amp + 1 h + 1 ppi
# 1 machine-amp + 1 free tRNA --> 1 machine + 1 amp + 1 charged tRNA
special_trna_subreactions = {}

initiation_subreactions = {
    'Translation_initiation_factor_InfA':
        {'enzymes': 'BSU01390-MONOMER',
         'stoich': {}},

    'Translation_initiation_factor_InfC':
        {'enzymes': 'BSU28870-MONOMER',
         'stoich': {}},

    'Translation_gtp_initiation_factor_InfB':
        {'enzymes': 'BSU16630-MONOMER',
         'stoich': {'gtp_c': -1,
                    'h2o_c': -1,
                    'h_c': 1,
                    'pi_c': 1,
                    'gdp_c': 1}},

    'fmet_addition_at_START':
        {'enzymes': ['BSU16630-MONOMER',
                     'BSU15730-MONOMER'],
         # iOL had h_c:1 for fmet addition but this is not mass balanced
         'stoich': {'10fthf_c': -1, 'thf_c': 1,
                    # 'h_c': 1,
                    'generic_tRNA_START_met__L_c': -1},
         'element_contribution': {'C': 1, 'O': 1}}
   }

elongation_subreactions = {'FusA_mono_elongation': {'enzymes': ['BSU01120-MONOMER'],
                                                    'stoich': {'gtp_c': -1,
                                                               'h2o_c': -1,
                                                               'h_c': 1,
                                                               'pi_c': 1,
                                                               'gdp_c': 1}},

                           'Tuf_gtp_regeneration': {'enzymes': ['BSU16500-MONOMER'],
                                                    'stoich': {}}}

# TODO Go through and double check elongation/termination ATP usage etc.

# The chaperones do not actively
# promote the folding of the substrate protein, but instead prevent ag-
# gregation of unfolded peptides. For a population of polypeptides,
# some fraction of the polypeptides released at the end of the cycle are in the
# native conformation. The remainder are rebound by DnaK or are di-
# verted to the chaperonin system (GroEL; see Fig. 4–31).

termination_subreactions = {'PrfA_mono_mediated_termination':
                            {'enzymes': ['BSU37010-MONOMER'],
                             'stoich': {}},

                            'PrfB_mono_mediated_termination':
                            {'enzymes': ['BSU35290-MONOMER'],
                             'stoich': {}},

                            'generic_RF_mediated_termination':
                            {'enzymes': ['generic_RF'],
                             'stoich': {}},

                            'N_terminal_methionine_cleavage':
                            {'enzymes': ['BSU01380-MONOMER'],
                             'stoich': {'h2o_c': -1,
                                        'met__L_c': 1, 'h_c': 1},
                             'element_contribution': {'H': -10, 'O': -1,
                                                      'C': -5, 'N': -1,
                                                      'S': -1}},

                            'peptide_deformylase_processing':
                            {'enzymes': ['BSU15720-MONOMER'],
                             'stoich': {'h2o_c': -1,
                                        'for_c': 1},
                             'element_contribution':
                                 {'H': 1, 'O': -1, 'C': -1}},

                            # This is a GTPS
                            'peptide_chain_release':
                            {'enzymes': ['BSU35290-MONOMER'],
                             'stoich': {'gtp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'pi_c': 1,
                                        'gdp_c': 1}},

                            'ribosome_recycler':
                            {'enzymes': ['BSU16520-MONOMER'],
                             'stoich': {}},

                            'GroEL_dependent_folding':
                            {'enzymes': ['CPLX8J2-24'],
                             'stoich': {'atp_c': -7,
                                        'h2o_c': -7,
                                        'h_c': 7,
                                        'adp_c': 7,
                                        'pi_c': 7}},

                            # DnaK is correct
                            'DnaK_dependent_folding':
                            {'enzymes': ['BSU25470-MONOMER', 'BSU25460-MONOMER',
                                         'BSU25480-MONOMER'],
                             'stoich': {'atp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'adp_c': 1,
                                        'pi_c': 1}}
                            }

# Subreaction for translation termination
translation_stop_dict = {'UAG': 'BSU37010-MONOMER',
                         'UGA': 'BSU35290-MONOMER',
                         'UAA': 'generic_RF'}

translation_start_codons = {"AUG", "GUG", "UUG", "AUU", "CUG"}

# Dictionary of frame shift mutations
frameshift_dict = {}

peptide_processing_subreactions = {"peptide_deformylase_processing",
                                   "peptide_chain_release",
                                   "ribosome_recycler"}


def add_translation_subreactions_to_model(me_model):
    # add general subreactions
    for rxn, info in elongation_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']

    # add subreactions associated with termination and postprocessing
    for rxn, info in termination_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})

    # add subreactions associated with translation initiation
    for rxn, info in initiation_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})


def add_charged_trna_subreactions(me_model):
    # create subreaction for each codon. this will be used to model
    # the addition of charged tRNAs to the elongating peptide
    for codon in dogma.codon_table:
        if dogma.codon_table[codon] == '*':
            stop_codon = codon.replace('T', 'U')
            stop_enzyme = translation_stop_dict.get(stop_codon)
            me_model.add_metabolites([Complex(stop_enzyme)])

            subreaction_data = SubreactionData(
                stop_codon + '_' + stop_enzyme + '_mediated_termination',
                me_model)
            subreaction_data.enzyme = stop_enzyme
            subreaction_data.stoichiometry = {}
        else:
            full_aa = dogma.amino_acids[dogma.codon_table[codon]]
            amino_acid = full_aa.split('_')[0]
            subreaction_data = SubreactionData(
                amino_acid + '_addition_at_' + codon.replace('T', 'U'),
                me_model)
            trna = 'generic_tRNA_' + codon.replace('T', 'U') + '_' + full_aa
            subreaction_data.enzyme = 'generic_Tuf'  # Default AA loader enzyme

            # Accounts for GTP hydrolyzed by EF-TU and the ATP hydrolysis to
            # AMP required to add the amino acid to the tRNA
            subreaction_data.stoichiometry = {'gtp_c': -1, 'h2o_c': -2,
                                              'gdp_c': 1, 'h_c': 2, 'pi_c': 1,
                                              'ppi_c': 1, 'amp_c': 1,
                                              'atp_c': -1,
                                              trna: -1}

    # Add subreactions for start codon and selenocysteine
    for rxn, info in special_trna_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})


# N terminal methionine cleaved
methionine_cleaved = ['BSU17410','BSU10790','BSU32150','BSU17460','BSU06480',
                      'BSU39920','BSU06180','BSU16100','BSU01390','BSU28230',
                      'BSU16150','BSU31390','BSU15990','BSU01440','BSU30540',
                      'BSU01290','BSU18000','BSU06850','BSU33540','BSU13900',
                      'BSU03520','BSU28310','BSU02890','BSU25020','BSU06150',
                      'BSU01050','BSU00730','BSU19530','BSU01410','BSU23860',
                      'BSU34790','BSU33940','BSU01100','BSU04730','BSU30190',
                      'BSU23040','BSU03130','BSU25410','BSU32710','BSU32890',
                      'BSU28440','BSU29120','BSU28500','BSU01340','BSU01020',
                      'BSU16250','BSU23470','BSU30650','BSU38550','BSU33910',
                      'BSU29660','BSU28870','BSU07000','BSU06030','BSU16170',
                      'BSU16690','BSU01310','BSU01150','BSU01700','BSU00510',
                      'BSU01200','BSU01040','BSU30120','BSU14580','BSU33400',
                      'BSU04190','BSU19550','BSU12290','BSU19240','BSU01120',
                      'BSU38140','BSU21870','BSU37660','BSU00110','BSU18030',
                      'BSU25480','BSU40030','BSU13180','BSU01780','BSU35000',
                      'BSU28430','BSU08820','BSU31350','BSU01250','BSU16500',
                      'BSU04400','BSU07830','BSU36830','BSU27320','BSU13470',
                      'BSU08070','BSU00480','BSU13160','BSU33900','BSU23820',
                      'BSU01220','BSU14600','BSU18239','BSU15080','BSU29490',
                      'BSU40890','BSU28860','BSU29080','BSU16680','BSU39340',
                      'BSU16520','BSU01300','BSU01190','BSU01030','BSU18360',
                      'BSU02900','BSU25550','BSU16490','BSU16330','BSU01110',
                      'BSU16180','BSU01500','BSU07850','BSU08550','BSU00520',
                      'BSU26660','BSU37540','BSU14590','BSU25470']
## GroEL dictionary from Endo & Kurusu, 2007
folding_dict = {
    'GroEL_dependent_folding': ['BSU23260','BSU10200','BSU28670','BSU14590',
                                'BSU23510','BSU05700','BSU15930','BSU27930',
                                'BSU25470','BSU28230','BSU14600','BSU14610',
                                'BSU36810','BSU36830','BSU01130','BSU33900',
                                'BSU29130','BSU16940','BSU31930','BSU33940',
                                'BSU35360','BSU29120','BSU16500','BSU25530',
                                'BSU37120','BSU39680','BSU11690','BSU40090'],
    'DnaK_dependent_folding': [
                              ]}

# Codons are not unique to a tRNA
trna_to_codon = {'BSU_tRNA_1': ['UUU', 'UUC'],
                 'BSU_tRNA_10': ['AUG','START'],
                 'BSU_tRNA_11': ['GAA', 'GAG'],
                 'BSU_tRNA_12': ['GUU', 'GUC', 'GUG', 'GUA'],
                 'BSU_tRNA_13': ['ACC', 'ACA', 'ACG', 'ACU'],
                 'BSU_tRNA_14': ['AAG', 'AAA'],
                 'BSU_tRNA_15': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_16': ['GGU', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_17': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_18': ['AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGU'],
                 'BSU_tRNA_19': ['CCG', 'CCA', 'CCU', 'CCC'],
                 'BSU_tRNA_2': ['GAU', 'GAC'],
                 'BSU_tRNA_20': ['GCA', 'GCG', 'GCC', 'GCU'],
                 'BSU_tRNA_21': ['AUG','START'],
                 'BSU_tRNA_22': ['GAU', 'GAC'],
                 'BSU_tRNA_23': ['AAC', 'AAU'],
                 'BSU_tRNA_24': ['ACC', 'ACA', 'ACG', 'ACU'],
                 'BSU_tRNA_25': ['GGU', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_26': ['AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGU'],
                 'BSU_tRNA_27': ['CCG', 'CCA', 'CCU', 'CCC'],
                 'BSU_tRNA_28': ['GCA', 'GCG', 'GCC', 'GCU'],
                 'BSU_tRNA_29': ['AAC', 'AAU'],
                 'BSU_tRNA_3': ['GAA', 'GAG'],
                 'BSU_tRNA_30': ['AGC', 'AGU', 'UCU', 'UCG', 'UCC', 'UCA'],
                 'BSU_tRNA_31': ['GAA', 'GAG'],
                 'BSU_tRNA_32': ['GUU', 'GUC', 'GUG', 'GUA'],
                 'BSU_tRNA_33': ['AUG','START'],
                 'BSU_tRNA_34': ['GAU', 'GAC'],
                 'BSU_tRNA_35': ['UUU', 'UUC'],
                 'BSU_tRNA_36': ['ACC', 'ACA', 'ACG', 'ACU'],
                 'BSU_tRNA_37': ['UAU', 'UAC'],
                 'BSU_tRNA_38': ['UGG'],
                 'BSU_tRNA_39': ['CAC', 'CAU'],
                 'BSU_tRNA_4': ['AAG', 'AAA'],
                 'BSU_tRNA_40': ['CAG', 'CAA'],
                 'BSU_tRNA_41': ['GGU', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_42': ['UGU', 'UGC'],
                 'BSU_tRNA_43': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_44': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_45': ['AAC', 'AAU'],
                 'BSU_tRNA_46': ['AGC', 'AGU', 'UCU', 'UCG', 'UCC', 'UCA'],
                 'BSU_tRNA_47': ['GAA', 'GAG'],
                 'BSU_tRNA_48': ['CAG', 'CAA'],
                 'BSU_tRNA_49': ['AAG', 'AAA'],
                 'BSU_tRNA_5': ['AUC', 'AUA', 'AUU'],
                 'BSU_tRNA_50': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_51': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_52': ['GUU', 'GUC', 'GUG', 'GUA'],
                 'BSU_tRNA_53': ['ACC', 'ACA', 'ACG', 'ACU'],
                 'BSU_tRNA_54': ['AAG', 'AAA'],
                 'BSU_tRNA_55': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_56': ['GGU', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_57': ['CUU', 'CUG', 'CUA', 'CUC', 'UUG', 'UUA'],
                 'BSU_tRNA_58': ['AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGU'],
                 'BSU_tRNA_59': ['CCG', 'CCA', 'CCU', 'CCC'],
                 'BSU_tRNA_6': ['GCA', 'GCG', 'GCC', 'GCU'],
                 'BSU_tRNA_60': ['GCA', 'GCG', 'GCC', 'GCU'],
                 'BSU_tRNA_61': ['AUG','START'],
                 'BSU_tRNA_62': ['AUG','START'],
                 'BSU_tRNA_63': ['AGC', 'AGU', 'UCU', 'UCG', 'UCC', 'UCA'],
                 'BSU_tRNA_64': ['AUG','START'],
                 'BSU_tRNA_65': ['GAU', 'GAC'],
                 'BSU_tRNA_66': ['UUU', 'UUC'],
                 'BSU_tRNA_67': ['CAC', 'CAU'],
                 'BSU_tRNA_68': ['GGU', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_69': ['AUC', 'AUA', 'AUU'],
                 'BSU_tRNA_7': ['AGC', 'AGU', 'UCU', 'UCG', 'UCC', 'UCA'],
                 'BSU_tRNA_70': ['AAC', 'AAU'],
                 'BSU_tRNA_71': ['AGC', 'AGU', 'UCU', 'UCG', 'UCC', 'UCA'],
                 'BSU_tRNA_72': ['GAA', 'GAG'],
                 'BSU_tRNA_73': ['CAG', 'CAA'],
                 'BSU_tRNA_74': ['CAG', 'CAA'],
                 'BSU_tRNA_75': ['GAA', 'GAG'],
                 'BSU_tRNA_76': ['ACC', 'ACA', 'ACG', 'ACU'],
                 'BSU_tRNA_77': ['UAU', 'UAC'],
                 'BSU_tRNA_78': ['GUU', 'GUC', 'GUG', 'GUA'],
                 'BSU_tRNA_79': ['AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGU'],
                 'BSU_tRNA_8': ['AUC', 'AUA', 'AUU'],
                 'BSU_tRNA_80': ['GGU', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_81': ['AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGU'],
                 'BSU_tRNA_82': ['GGU', 'GGG', 'GGA', 'GGC'],
                 'BSU_tRNA_83': ['GUU', 'GUC', 'GUG', 'GUA'],
                 'BSU_tRNA_84': ['AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGU'],
                 'BSU_tRNA_85': ['AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGU'],
                 'BSU_tRNA_86': ['GCA', 'GCG', 'GCC', 'GCU'],
                 'BSU_tRNA_9': ['GCA', 'GCG', 'GCC', 'GCU']}

from __future__ import division, absolute_import, print_function

from six import iteritems

from cobrame import ComplexData, SubreactionData

transcription_subreactions = {
    'Transcription_normal_rho_independent':
        {'enzymes': ['BSU00550-MONOMER', 'BSU16600-MONOMER', 'BSU01010-MONOMER',
                     'BSU27320-MONOMER','CPLX8J2-30'],
         'stoich': {}},
    'Transcription_normal_rho_dependent':
        {'enzymes': ['BSU00550-MONOMER', 'BSU16600-MONOMER', 'BSU01010-MONOMER',
                     'BSU27320-MONOMER','CPLX8J2-30', 'BSU37080-MONOMER'],
         'stoich': {'atp_c': -3,
                    'h2o_c': -3,
                    'adp_c': 3,
                    'pi_c': 3,
                    'h_c': 3}},
    'Transcription_stable_rho_independent':
        {'enzymes': ['BSU00550-MONOMER', 'BSU16600-MONOMER', 'BSU01010-MONOMER',
                     'BSU27320-MONOMER','CPLX8J2-30'],
         'stoich': {}},
    'Transcription_stable_rho_dependent':
        {'enzymes': ['BSU00550-MONOMER', 'BSU16600-MONOMER', 'BSU01010-MONOMER',
                     'BSU27320-MONOMER','CPLX8J2-30', 'BSU37080-MONOMER'],
         'stoich': {'atp_c': -3,
                    'h2o_c': -3,
                    'adp_c': 3,
                    'pi_c': 3,
                    'h_c': 3}}
}

sigma_factor_complex_to_rna_polymerase_dict = {
      'BSU25200-MONOMER': 'CPLX8J2-52', # sigA (rpoD)
      'BSU15320-MONOMER': 'CPLX8J2-61', # sigE
      'BSU15330-MONOMER': 'CPLX8J2-60', # sigG
      'BSU23450-MONOMER': 'CPLX8J2-58', # sigF
      'BSU34200-MONOMER': 'CPLX8J2-57', # sigL
      'BSU04730-MONOMER': 'CPLX8J2-56', # sigB
      'BSU00980-MONOMER': 'CPLX8J2-55', # sigH
      'BSU16470-MONOMER': 'CPLX8J2-54', # sigD
      'BSU13450-MONOMER': 'CPLX8J2-53', # sigI
      'BSU14730-MONOMER': 'CPLX8J2-51', # sigECF
      'BSU01730-MONOMER': 'CPLX8J2-50', # sigW
      'BSU23100-MONOMER': 'CPLX8J2-49', # sigX
      'BSU26840-MONOMER': 'CPLX8J2-48', # sigZ
      'BSU27120-MONOMER': 'CPLX8J2-47', # sigV
      'BSU38700-MONOMER': 'CPLX8J2-46', # sigY
      'BSU09520-MONOMER': 'CPLX8J2-45', # sigM
      'MONOMER8J2-6': 'CPLX8J2-60', # sigK
      'CPLX8J2-36': 'CPLX8J2-30' # 2-subunit
      }

rna_polymerase_id_by_sigma_factor = {} # Complexes are already formed

rna_polymerases = list(rna_polymerase_id_by_sigma_factor.keys())

# Degradosome composition from Lehnik-Habrink, 2010.
rna_degradosome = {'BSU16960-MONOMER': 1, 'CPLX8J2-39': 1,'BSU16690-MONOMER': 1,
                   'BSU33900-MONOMER': 1, 'BSU29190-MONOMER': 1}

## BSU41050-MONOMER is set instead of CPLX8J2-62, since it requires RNA thats no longer available after pruning.
excision_machinery = {
    'rRNA_containing': ['BSU14530-MONOMER','BSU15930-MONOMER','BSU00410-MONOMER'],
    'monocistronic': ['BSU41050-MONOMER', 'BSU23840-MONOMER',
                      'generic_RNase'],
    'polycistronic_wout_rRNA': ['BSU41050-MONOMER', 'BSU23840-MONOMER',
                                'generic_RNase']}


def add_rna_polymerase_complexes(me_model, verbose=True):

    for cplx, components in iteritems(rna_polymerase_id_by_sigma_factor):
        rnap_complex = ComplexData(cplx, me_model)
        rnap_components = rnap_complex.stoichiometry
        sigma_factor = components['sigma_factor']
        polymerase = components['polymerase']

        rnap_components[sigma_factor] = 1
        rnap_components[polymerase] = 1

        rnap_complex.create_complex_formation(verbose=verbose)


def add_rna_splicing(me_model):

    # Ecoli has three alternatie mechanisms for splicing RNA, depending
    # on what RNA types the TU contains
    excision_types = ['rRNA_containing', 'monocistronic',
                      'polycistronic_wout_rRNA']

    for excision_type in excision_types:
        complex_data = ComplexData(excision_type + "_excision_machinery",
                                   me_model)

        for machine in excision_machinery[excision_type]:
            complex_data.stoichiometry[machine] = 1

        complex_data.create_complex_formation()
        modification = SubreactionData(excision_type + "_excision", me_model)
        modification.enzyme = complex_data.id

    # Loop through transcription reactions and add appropriate splicing
    # machinery based on RNA types and number of splices required
    for t in me_model.transcription_data:
        n_excised = sum(t.excised_bases.values())
        n_cuts = len(t.RNA_products) * 2
        if n_excised == 0 or n_cuts == 0:
            continue
        rna_types = list(t.RNA_types)
        n_trna = rna_types.count("tRNA")

        if "rRNA" in set(rna_types):
            t.subreactions["rRNA_containing_excision"] = n_cuts
        elif n_trna == 1:
            t.subreactions["monocistronic_excision"] = n_cuts
        elif n_trna > 1:
            t.subreactions["polycistronic_wout_rRNA_excision"] = n_cuts
        else:  # only applies to rnpB
            t.subreactions["monocistronic_excision"] = n_cuts

        # The non functional RNA segments need degraded back to nucleotides
        # TODO check if RNA_degradation requirement is per nucleotide
        t.subreactions["RNA_degradation_machine"] = n_cuts
        t.subreactions["RNA_degradation_atp_requirement"] = n_excised

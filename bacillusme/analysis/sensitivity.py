from __future__ import print_function, division, absolute_import

import sys

import qminospy
from qminospy.me2 import ME_NLP

# python imports
from copy import copy
import re
from os.path import join
from collections import defaultdict
import pickle

# third party imports
import pandas
import cobra
from tqdm import tqdm
import numpy as np
import scipy

# COBRAme
import cobrame
from cobrame.util import building, mu, me_model_interface
from cobrame.io.json import save_json_me_model, save_reduced_json_me_model

# ECOLIme
import bacillusme
from bacillusme import (transcription, translation, flat_files, generics, formulas, compartments)
from bacillusme.util.helper_functions import *
from cobrame.core.model import MEModel


# Sensitivity functions

def add_dummy_demand(me,met_id,flux=0,met_flux_dict=0,fraction=0.01):
    '''
    This function adds a dummy demand reactions to calculate sensitivity and cost
    '''
    met = me.metabolites.get_by_id(met_id)
    rxn_id = 'dummy_demand'
    # Replace dummy_demand if existent
    if hasattr(me.reactions,rxn_id):
        me.reactions.get_by_id(rxn_id).remove_from_model()
    
    if isinstance(me,MEModel):
        rxn = cobrame.MEReaction(rxn_id)
    else:
        rxn = cobra.Reaction(rxn_id)
    me.add_reaction(rxn)
    
    if isinstance(met,cobrame.Metabolite) or not isinstance(me,MEModel):
        rxn.reaction = met_id + ' ->'
        mass_correction = 0
    elif isinstance(met,cobrame.Complex):
        mass_correction = 'protein_biomass'
        weight = met.formula_weight/1000 # g/mmol
        rxn.reaction = met_id + ' + ' + str(weight) + ' ' + mass_correction + ' ->'
    else:
        raise ValueError('Metabolite type not supported')
    
    if met_flux_dict and not flux:
        flux = met_flux_dict[met_id]
    
    if not flux and mass_correction:
        # fraction 1% of total group production. This is done if production of a metabolite causes the
        # production of biomass, e.g. protein_biomass, lipid_biomass, etc.
        group_production = me.solution.x_dict[mass_correction+'_to_biomass']
        flux = fraction * group_production / weight
    elif not flux:
        flux = 0.001
        
    rxn.lower_bound = flux
    rxn.upper_bound = flux
    
    
def single_change(me,met_id,single_change_function=False,met_flux_dict=0,growth_key='mu'):
    if single_change_function=='transporter': # close_transporter function
        close_transporter(me,met_id)
    elif single_change_function == 'overexpress':
        overexpress_transporter(me,met_id)
    elif single_change_function == 'group_knockout':
        group_knockout(me,met_id)
    elif single_change_function == 'feed_metabolite':
        feed_metabolite(me,met_id,met_flux_dict=met_flux_dict)
    elif single_change_function == 'gene_knockout':
        single_gene_knockout(me, met_id,growth_key=growth_key)
    elif single_change_function == 'knockout_protein_sink':
        knockout_protein_sink(me,met_id,growth_key=growth_key)
    else: # Just normal sensitivity/cost calculation
        add_dummy_demand(me,met_id,met_flux_dict=met_flux_dict)
    
def single_flux_response(me,met_id,mu_fix=False,precision=1e-6,\
                         single_change_function=False, growth_key='mu', met_flux_dict=0):
    '''
    This function calculates flux response of a single metabolite. Response of growth and energy are
    sensitixvity and cost, respectively.
    '''

    if isinstance(met_id,list):
        for met in met_id:
            single_change(me,met,single_change_function=single_change_function,growth_key=growth_key)
        met_id = met
    else:
        single_change(me,met_id,single_change_function=single_change_function,met_flux_dict=met_flux_dict,growth_key=growth_key)
        
    solve_model(me,precision = precision,mu_fix=mu_fix,growth_key=growth_key)
    
    if me.solution and me.solution.f > 0.:
        flux_dict = me.solution.x_dict
    else:
        flux_dict = {r.id:0. for r in me.reactions}
    return met_id, flux_dict

def solve_model(model,precision = 1e-6,mu_fix=False,growth_key='mu'):
    if isinstance(model,cobrame.core.model.MEModel):
        solve_me_model(model, 0.2, min_mu = .05, using_soplex=False,precision = precision,mu_fix=mu_fix,
                       verbosity=0,growth_key=growth_key)
    else:
        model.optimize()
        

def all_flux_responses(me,met_ids,mu_fix=False,solution=0,NP=1,precision=1e-6,
                       single_change_function=False,growth_key = 'mu',sequential=False,met_flux_dict=0):
    '''
    This function calculates flux responses of several metabolites. Response of growth and energy are
    sensitivity and cost, respectively.
    '''
    import copy
    # Correct number of threads if necessary
    NP = min([NP,len(met_ids)])
    if not solution:
        solve_model(me,precision = precision,mu_fix=mu_fix)
    flux_dict = dict()
    flux_dict['base'] = me.solution.x_dict

    obj = [rxn.id for rxn in me.objective][0]
    
    
    if NP == 1:
        for met_id in tqdm(met_ids,leave=True,position=0):
            og_me = copy.deepcopy(me)
            single_flux_response(og_me,met_id,single_change_function=single_change_function,\
                                 mu_fix=mu_fix,precision=precision,growth_key = growth_key,met_flux_dict=met_flux_dict)
            flux_dict[met_id] = og_me.solution.x_dict
    else:
        import multiprocessing as mp
        
        pool = mp.Pool(NP)
        pbar = tqdm(total=len(met_ids),position=0,leave=True)
        pbar.set_description('{} response ({} threads)'.format(obj,NP))
        def collect_result(result):
            pbar.update(1)
            flux_dict[result[0]] = result[1]

        for met_id in met_ids:
            if sequential:
                met_idx = met_ids.index(met_id)
                met_arg = met_ids[:met_idx+1] # All mets until met
            else:
                met_arg = met_id
            args = (me,met_arg)
            kwds = {'single_change_function':single_change_function,\
                'mu_fix':mu_fix,'precision':precision,'growth_key':growth_key,
                'met_flux_dict':met_flux_dict}
            pool.apply_async(single_flux_response,args,kwds,
                callback=collect_result)
        pool.close()
        pool.join()
    flux_results_df = pd.DataFrame.from_dict(flux_dict)
    return flux_results_df

def process_flux_responses(me,flux_results_df,base_rxn):
    base_flux = flux_results_df.loc[base_rxn]['base']
    met_ids = flux_results_df.columns.values
    processed = {}
    for met_id in met_ids:
        if met_id == 'base':
            continue
        processed[met_id] = {}
        new_flux = flux_results_df.loc[base_rxn][met_id]
        met = me.metabolites.get_by_id(met_id)
        met_flux = flux_results_df.loc['dummy_demand'][met_id]
        response = (base_flux-new_flux)/met_flux
        try:
            MW = met.formula_weight/1000
        except:
            MW = 0
        processed[met_id]['molecular_weight'] = MW
        processed[met_id]['sensitivity'] = response
    processed_results_df = pd.DataFrame.from_dict(processed).T

    return processed_results_df

def sensitivity(me,met_ids,NP=1,solution=0,biomass_dilution='biomass_dilution',\
		precision=1e-6,met_flux_dict=False):
	'''
	This function calculates sensitivity of all metabolites. 
	Parallel processing is activated when NP>1.
	'''

	flux_results_df = all_flux_responses(me,met_ids,mu_fix=False,\
			solution=solution,NP=NP,precision=1e-6,met_flux_dict=met_flux_dict)

	# Process results to get sensitivity
	sensitivity_df = process_flux_responses(me,flux_results_df,biomass_dilution)

	return sensitivity_df,flux_results_df

def biosynthetic_cost(me,met_ids,cost_rxn='ATPM',\
		biomass_dilution='biomass_dilution', NP=1,solution=0,precision=1e-6):
	'''
	This function calculates sensitivity of all metabolites. 
	Parallel processing is activated when NP>1.
	'''
	if not solution:
		solve_model(me,precision = precision,mu_fix=False,growth_key='mu')
	base_flux = me.solution.x_dict[cost_rxn]
	base_mu = me.solution.f

	# Growth rate-fixed solution.
	me.objective = cost_rxn
	mu_fix = 0.9*base_mu
	if not isinstance(me,MEModel):
		me.reactions.get_by_id(biomass_dilution).bounds = (mu_fix,mu_fix)

	flux_results_df = all_flux_responses(me,met_ids,NP=NP,\
			mu_fix=mu_fix,solution=0,precision=1e-6)

	# Process results to get cost
	cost_df = process_flux_responses(me,flux_results_df,cost_rxn)

	return cost_df,flux_results_df


def close_transporter(me,transport_id):
    r = me.reactions.get_by_id(transport_id)
    r.lower_bound = 0
    r.upper_bound = 0
        
def overexpress_transporter(me,transport_id):
    r = me.reactions.get_by_id(transport_id)
    base_flux = me.solution.x_dict[r.id]
        
    r.lower_bound = base_flux*2
    r.upper_bound = base_flux*2
        
def group_knockout(me,met_id):
    transport_reactions = get_transport_reactions(me,met_id,comps=['c','s']) \
                + get_transport_reactions(me,met_id,comps=['s','c'])
    
    for r in transport_reactions:
        if 'SPONT' not in r.id:
            r.lower_bound = 0
            r.upper_bound = 0    
            
def partial_knockout(me,met_id,frac=0.5):
    transport_reactions = get_transport_reactions(me,met_id,comps=['c','s']) \
                + get_transport_reactions(me,met_id,comps=['s','c'])
    
    for r in transport_reactions:
        if 'SPONT' not in r.id:
            r.upper_bound = frac * me.solution.x_dict[r]
            
def feed_metabolite(me,met_id,met_flux_dict=0):
    fed_met = re.split('_[0-9]*$',met_id)[0] +'_e'
    rxn_id = 'EX_{}'.format(fed_met)

    if not hasattr(me.reactions,rxn_id):
        return
    
    met = me.metabolites.get_by_id(fed_met)

    rxn = me.reactions.get_by_id(rxn_id)
    rxn.upper_bound = 0
    
    if met_flux_dict:
        rxn.lower_bound = met_flux_dict[met_id]
        rxn.upper_bound = met_flux_dict[met_id]
    else:
        rxn.lower_bound = -1000
    
#     print(rxn.reaction,rxn.lower_bound,rxn.upper_bound)

def knockout_protein_sink(me,met_id,growth_key='mu'):
    single_gene_knockout(me, met_id,precision=1e-6,growth_key=growth_key)
    me.reactions.get_by_id('SK_protein_{}'.format(met_id)).bounds=(0,0)

def transporter_knockout(me,transport_ids,NP=1,solution=0,precision=1e-6,growth_key='mu',\
                        biomass_dilution='biomass_dilution',single_change_function='transporter',sequential=False):
    '''
    This function calculates the response of shutting down
    transporters.
    ''' 
    
    print('Chosen change function: {}      Sequential = {}'.format(single_change_function,str(sequential)))
        
    flux_results_df = all_flux_responses(me,transport_ids,mu_fix=False,\
            solution=solution,NP=NP,precision=1e-6,\
            single_change_function=single_change_function,growth_key=growth_key,sequential=sequential)

    return flux_results_df

def single_gene_knockout(model, gene_id,precision=1e-6,growth_key='mu'):
    temp_model = model
    if not isinstance(model,cobrame.MEModel):
        temp_gene = temp_model.genes.get_by_id(gene_id)
        reactions = temp_gene.reactions
        for reaction in reactions:
            rule = reaction.gene_reaction_rule
            individual_rules = rule.split(' or ')
            if len(individual_rules) == 1:
                reaction.lower_bound = 0.
                reaction.upper_bound = 0.
        temp_model.optimize()
    elif isinstance(model, cobrame.MEModel):
        temp_model.reactions.get_by_id('translation_' + gene_id).bounds = (0,0)
        solve_me_model(temp_model, 1., min_mu = .1,
                       precision=precision, using_soplex=False,
                       verbosity=0,growth_key=growth_key)

    if model.solution and model.solution.status == 'optimal':
        flux_dict = model.solution.x_dict
    else:
        flux_dict = {r.id:0. for r in model.reactions}

    return gene_id, flux_dict
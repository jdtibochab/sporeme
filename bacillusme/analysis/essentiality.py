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

def single_gene_knockout(model, gene_id, model_type,precision=1e-6):
	temp_model = copy(model)
	if model_type == 'm':
		temp_model = model.copy()
		temp_gene = temp_model.genes.get_by_id(gene_id)
		reactions = temp_gene.reactions
		for reaction in reactions:
			rule = reaction.gene_reaction_rule
			individual_rules = rule.split(' or ')
			if len(individual_rules) == 1:
				reaction.lower_bound = 0.
				reaction.upper_bound = 0.
		temp_model.optimize()
	elif model_type == 'me':
		temp_model = copy(model)
		protein = temp_model.metabolites.get_by_id('protein_' + gene_id)
		reactions = protein.reactions
		for reaction in reactions:
			reaction.lower_bound = 0.
			reaction.upper_bound = 0.
		solve_me_model(temp_model, 1., min_mu = .1, precision=precision, using_soplex=False,verbosity=0)

	if temp_model.solution and temp_model.solution.status == 'optimal':
		flux_dict = temp_model.solution.x_dict
	else:
		flux_dict = {r.id:0. for r in temp_model.reactions}
        
	return gene_id, flux_dict

def gene_knockout_responses(model,genes, model_type, NP = 1, precision=1e-6):
    import numpy as np
    flux_dict = {}

    ## Calculation
    NP = min([NP,20,len(genes)])
    if NP == 1:
        from copy import copy
        for gene in tqdm(genes):
            _,single_flux_dict = single_gene_knockout(model, gene, model_type)
            flux_dict[gene] = single_flux_dict
    else:
        
        pbar = tqdm(total=len(genes))
        pbar.set_description('{} threads'.format(NP))
        def collect_result(result):
            pbar.update(1)
            flux_dict[result[0]] = result[1]
        import multiprocessing as mp
        # Initiate pool
        pool = mp.Pool(NP)
        # Calculation
        for gene in genes:
            args = (model,gene,model_type)
            kwds = {'precision':precision}
            pool.apply_async(single_gene_knockout,args,kwds,
                callback=collect_result)
        # Close
        pool.close()
        pool.join()
    return flux_dict

def gene_essentiality(model,genes, model_type,  threshold = 0.01, NP = 1,
                        initial_f = 0,precision=1e-6, flux_responses=None):
    ## Initialization
    if isinstance(flux_responses,type(None)):
        if not initial_f:
            print('Solving model for base solution')
            if model_type == 'm':
                model.optimize()
            elif model_type =='me':
                solve_me_model(model, max_mu=1., min_mu = .1, precision=precision, using_soplex=False, verbosity=0)
            if model.solution.status != 'optimal':
                raise Exception("Model not feasible")
            flux_dict = gene_knockout_responses(model,genes, model_type, NP=NP , precision=precision)
            flux_dict['base'] = model.solution.x_dict
            flux_responses = pd.DataFrame.from_dict(flux_dict)
            flux_responses[abs(flux_responses<1e-15)] = 0.

    growth_rxn = [rxn.id for rxn in model.objective][0]
    growths = flux_responses.loc[growth_rxn]
    response = (growths-growths['base'])/growths['base']

    essentiality = pd.Series(index=response.index)
    for gene,change in response.iteritems():
        if abs(change) < threshold:
            essentiality[gene] = '0' # No change
        elif change < threshold-1:
            essentiality[gene] = 'e' # Essential
        elif change < -threshold:
            essentiality[gene] = '-' # Growth reduced
        elif change > threshold:
            essentiality[gene] = '+' # Growth increased
        else:
            essentiality[gene] = 'e' # Model broke = essential

    frame = {'response':response,'essentiality':essentiality}
    essentiality_df = pd.DataFrame(frame).drop('base')
    print('Done')
    return essentiality_df.replace(to_replace=0.0, value='0'),flux_responses
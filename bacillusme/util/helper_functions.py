import cobrame
from cobrame.core.model import MEModel
from tqdm import tqdm
import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np

def get_base_complex_data(model, complex_id):
    """If a complex is modified in a metabolic reaction it will not
    have a formation reaction associated with it. This function returns
    the complex data of the "base" complex, which will have the subunit
    stoichiometry of that complex"""

    # First try unmodified complex id
    try_1 = complex_id.split('_')[0]
    if try_1 in model.process_data:
        return model.process_data.get_by_id(try_1)

    try_2 = complex_id.split('_')[0] + '_'
    count = 0
    for i in model.process_data.query(try_2):
        if isinstance(i, cobrame.ComplexData):
            count += 1
            data = i
    if count == 0:
        raise UserWarning('No base complex found for %s' % complex_id)
    if count > 1:
        raise UserWarning('More than one possible base complex found for %s' %
                          complex_id)
    return data

def get_identical_reactions(ref,rxn):
    candidate_rxns = []
    metabolites = rxn.metabolites
    met_ids = []
    # Get reaction metabolite IDs
    for metabolite in metabolites:
        met_ids.append(metabolite.id)
    # Look for identical reactions in reference model
    for ref_rxn in ref.reactions:
        ref_metabolites = ref_rxn.metabolites
        ref_met_ids = []
        if len(metabolites) == len(ref_metabolites):
            for ref_metabolite in ref_metabolites:
                ref_met_ids.append(ref_metabolite.id)
            if len(list(set(ref_met_ids) & set(met_ids))) == len(metabolites):
                candidate_rxns.append(ref_rxn)

    return candidate_rxns

def get_gene_info(gb_file,info,ID,element_types):
    output = None
    for feature in gb_file.features:
    # Skip if not a gene used in ME construction
        if feature.type not in element_types or 'pseudo' in feature.qualifiers:
            continue
        if feature.qualifiers["locus_tag"][0] == ID:
            output = feature.qualifiers[info]
    return output

    ############################### NEW ###############################

def test_metabolite_production(me,metabolites,muf = 0.):
    from qminospy.me2 import ME_NLP
    gap_mets = []

    if not muf and me.global_info['k_deg'] != 0:
        print ('Updating model with kdeg = 0 for mu = 0')
        me.global_info['k_deg'] = 0.
        me.update()

    for met_id in metabolites:

        r_id = 'DM_' + met_id
        r = cobrame.MEReaction(r_id)
        try:
            me.add_reaction(r)
            r.reaction = met_id + '->'
        except:
            print(me.reactions.get_by_id(r_id).id,' already in model')
        #print_reactions_of_met(me,met_id)
        me.objective = r_id

        me_nlp = ME_NLP(me, growth_key='mu')
        x,status,hs = me_nlp.solvelp(muf)

        f = me.solution.x_dict[r_id]

        if not status == 'optimal' or f < 0.01:
            gap_mets.append(met_id)
        print(met_id, status, f)
    return gap_mets

def identify_precursors(me,metabolite_id,only_direct_precursors = False,ignore_classes = None, force_classes = None):

    ### NOT WORKING YET

    import copy
    precursors = []
    formation_reactions = []
    metabolite = me.metabolites.get_by_id(metabolite_id)
    for rxn in me.reactions:
        if metabolite in rxn.products:
            formation_reactions.append(rxn)

    for rxn in formation_reactions:
        for reactant in rxn.reactants:
            if reactant.id not in precursors:
                precursors.append(reactant.id)
    direct_precursors = copy.copy(precursors)
    #print(precursors)
    if not only_direct_precursors:
        for precursor in precursors:
            for rxn in me.metabolites.get_by_id(precursor).reactions:
                products_of_rxn = [product.id for product in rxn.products]
                if precursor in products_of_rxn:
                    reactants = rxn.reactants
                    reactant_ids = [met.id for met in reactants]
                    if metabolite_id not in reactant_ids:            
                        for reactant in reactants:
                            if reactant.id not in precursors:
                                precursors.append(reactant.id)
                                #print(reactant.id)
    test_precursors = copy.copy(precursors)
    if ignore_classes:
        for precursor_id in test_precursors:
            precursor = me.metabolites.get_by_id(precursor_id)
            for ignore_class in ignore_classes:
                if isinstance(precursor,ignore_class):
                    precursors.remove(precursor_id)
                    break
    print(len(precursors))
    test_precursors = copy.copy(precursors)
    if force_classes:
        for precursor_id in test_precursors:
            precursor = me.metabolites.get_by_id(precursor_id)
            e = 1
            for force_class in force_classes:
                if isinstance(precursor,force_class):
                    e = 0
            if e:
                precursors.remove(precursor_id)	
    print(len(precursors))
    return precursors,direct_precursors

def get_reactions_of_met(me,met,s = 0, ignore_types = (),only_types = (), verbose = True):
    import copy
    met_stoich = 0
    if only_types:
        only_reaction_types = tuple([getattr(cobrame,i) for i in only_types])
    elif ignore_types:
        ignore_reaction_types = tuple([getattr(cobrame,i) for i in ignore_types])
    reactions = []

    if not hasattr(me.metabolites,met):
        return reactions
    for rxn in me.metabolites.get_by_id(met).reactions:
        if only_types and not isinstance(rxn, only_reaction_types):
            continue
        elif ignore_types and isinstance(rxn, ignore_reaction_types):
            continue
        reactants = [met.id for met in rxn.reactants]
        products = [met.id for met in rxn.products]

        try:
            pos = 1 if met in products else -1
            rev = 1 if rxn.lower_bound < 0 else 0
            fwd = 1 if rxn.upper_bound > 0 else 0
        except:
            if verbose:
                print(rxn.id, 'symbolic bounds')
            else:
                pass

        try:
            if not s:
                reactions.append(rxn)
                if verbose:
                    print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)

            elif s == pos*fwd or s == pos*rev:
                reactions.append(rxn)
                if verbose:
                    print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)

        except:
            if verbose:
                print(rxn.id, 'no reaction')
            else:
                pass
    return reactions

def add_exchange_reactions(me,metabolites):
    for met in metabolites:
        rxn_id = "EX_" + met
        try:
            r = cobrame.MEReaction(rxn_id)
            me.add_reaction(r)
            r.reaction = met + " <=> "
        except:
            r = me.reactions.get_by_id(rxn_id)
        r.lower_bound = -1000
        #print(r.id,r.lower_bound,r.upper_bound,r.reaction)
    return me

def brute_force_check(me,metabolites_to_add,objective_function = 'biomass_dilution',muf = 0.01, min_f = 0.01):

    me.objective = objective_function

    from qminospy.me2 import ME_NLP
    print('Added exchange reactions ')
    me = add_exchange_reactions(me,metabolites_to_add)
    print('Objective: ', objective_function, me.reactions.get_by_id(objective_function).reaction)

    me_nlp = ME_NLP(me, growth_key='mu')
    x,status,hs = me_nlp.solvelp(muf)
    initial_f = me.solution.x_dict[objective_function]
    print('Initial objective function value of ', initial_f, status)

    if not status =='optimal':
        return

    if initial_f < min_f:
        print('No production capacity of objective')

    print(me.solution.x_dict['formation_ribosome'])

    eliminate_mets = []
    for met_id in metabolites_to_add:
        ex_rxn_id = "EX_" + met_id
        ex_rxn_flux = me.solution.x_dict[ex_rxn_id]
        ex_rxn = me.reactions.get_by_id(ex_rxn_id)
        if ex_rxn_flux > 0:
            me.reactions.get_by_id(ex_rxn_id).lower_bound = 0
            me.reactions.get_by_id(ex_rxn_id).upper_bound = 1000
            print(ex_rxn_id, ex_rxn_flux, ex_rxn.reaction)
        elif ex_rxn_flux < 0:
            me.reactions.get_by_id(ex_rxn_id).lower_bound = -1000
            me.reactions.get_by_id(ex_rxn_id).upper_bound = 0
            print(ex_rxn_id, ex_rxn_flux, ex_rxn.reaction)
        elif ex_rxn_flux == 0:
            me.reactions.get_by_id(ex_rxn_id).lower_bound = 0
            me.reactions.get_by_id(ex_rxn_id).upper_bound = 0
            print(ex_rxn_id, ' carrying no flux ... eliminated')

            eliminate_mets.append(met_id)


    for el_met_id in eliminate_mets:
        el_rxn_id = 'EX_' + el_met_id
        metabolites_to_add.remove(el_met_id)

    print('Processing ', len(metabolites_to_add), ' metabolites')

    gap_mets = []
    for met_id in metabolites_to_add:
        ex_rxn_id =  "EX_" + met_id

        lb = me.reactions.get_by_id(ex_rxn_id).lower_bound
        ub = me.reactions.get_by_id(ex_rxn_id).lower_bound

        me.reactions.get_by_id(ex_rxn_id).lower_bound = 0
        me.reactions.get_by_id(ex_rxn_id).upper_bound = 0

        me_nlp = ME_NLP(me, growth_key='mu')
        x,status,hs = me_nlp.solvelp(muf)
        f = me.solution.x_dict[objective_function]

        el_bool = ''
        if not status == 'optimal' or f < min_f:
            me.reactions.get_by_id(ex_rxn_id).lower_bound = lb
            me.reactions.get_by_id(ex_rxn_id).lower_bound = ub
            gap_mets.append(met_id)
            el_bool = ' gap'

        print(met_id, status, f, el_bool, '... Gaps: ', len(gap_mets))
    return gap_mets

def solve_me_model(me, max_mu=1., precision=1e-6, min_mu=0, using_soplex=True,
                  compiled_expressions=None, verbosity = 2, mu_fix = False,
                  growth_key='mu'):
    from qminospy.me1 import ME_NLP1
    ## If fixed growth rate, solve as LP
    if mu_fix:
        me_nlp = ME_NLP1(me)
        me_nlp.solvelp(mu_fix)
    else:
        ## 
        if using_soplex:
            from cobrame.solve.algorithms import binary_search
            binary_search(me, min_mu=min_mu, max_mu=max_mu, debug=True, mu_accuracy=precision,
                compiled_expressions=compiled_expressions)
        else:
            # The object containing solveME methods--composite that uses a ME model object 
            me_nlp = ME_NLP1(me, growth_key=growth_key)
            # Use bisection for now (until the NLP formulation is worked out)
            muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=max_mu, verbosity=verbosity)
            try:
                me.solution.f = me.solution.x_dict['biomass_dilution']
            except:
                pass

def show_escher_map(me, solution=None):
    import escher
    view = escher.Builder("iJO1366.Central metabolism")
    view.reaction_data = me.get_metabolic_flux(solution=solution)
    return view

def open_all_exchange(me):
    for rxn in me.reactions:
        rxn_id = rxn.id
        if 'EX_' in rxn_id:
            rxn.upper_bound = 1000
            rxn.lower_bound = -1000
    return me

def is_same_reaction(rxn,ref_rxn,approximate=False):
    a = 2 if approximate else 0
    reactants = [met.id[0:len(met.id)-a] for met in rxn.reactants]
    ref_reactants = [met.id[0:len(met.id)-a] for met in ref_rxn.reactants]
    products = [met.id[0:len(met.id)-a] for met in rxn.products]
    ref_products = [met.id[0:len(met.id)-a] for met in ref_rxn.products]
    if set(reactants)==set(ref_reactants) and set(products)==set(ref_products):
        return 1
    else:
        return 0
    
def homogenize_metabolites(model,ref_model,approximate=False):
    
    # Warning: This function assumes metabolite IDs are conserved.
    met_dict = {}
    for met in model.metabolites:
        if approximate:
            met_dict[met.id] = [ref_met.id for ref_met in ref_model.metabolites\
                              if ref_met.id[0:len(ref_met.id)-2] == met.id[0:len(met.id)-2]]
        elif met.id in ref_model.metabolites:
            met_dict[met.id] = [met.id]
        else:
            met_dict[met.id] = []
    return met_dict
    

def homogenize_reactions(model,ref_model,approximate=False):
    print('Homogenizing metabolites {} against {}'.format(model.id,ref_model.id))
    met_dict = homogenize_metabolites(model,ref_model,approximate=approximate)
    all_ref_rxns = [rxn.id for rxn in ref_model.reactions]
    rxn_dict = dict()
    rxn_id_dict = {}
    
    print('Homogenizing reactions {} against {}'.format(model.id,ref_model.id))
    for rxn in tqdm(model.reactions):
        if rxn not in rxn_dict:
            rxn_dict[rxn] = []
            rxn_id_dict[rxn.id] = []
        base_met = rxn.reactants[0]
        for m in met_dict[base_met.id]:
            ref_met = ref_model.metabolites.get_by_id(m)
            for ref_rxn in ref_met.reactions:
                if is_same_reaction(rxn,ref_rxn,approximate=approximate):
                    rxn_dict[rxn].append(ref_rxn)
                    rxn_id_dict[rxn.id].append(ref_rxn.id)
    return rxn_dict,rxn_id_dict

def exchange_single_model(me, flux_dict = 0, solution=0):
    import pandas as pd

    complete_dict = {'id':[],'name':[],'reaction':[],'lower_bound':[],'upper_bound':[],'flux':[]}

    if solution:
        flux_dict = solution.x_dict
    elif not flux_dict:
        flux_dict = me.solution.x_dict

    for rxn in me.reactions:
        if 'EX_' in rxn.id:
            flux = flux_dict[rxn.id]

            if not flux:
                continue
            rxn_name = rxn.name
            reaction = rxn.reaction
            lb = rxn.lower_bound
            ub = rxn.upper_bound

            complete_dict['id'].append(rxn.id)
            complete_dict['name'].append(rxn_name)
            complete_dict['reaction'].append(reaction)
            complete_dict['lower_bound'].append(lb)
            complete_dict['upper_bound'].append(ub)
            complete_dict['flux'].append(flux)


    df = pd.DataFrame(complete_dict).set_index('id')
    return df

def get_metabolites_from_pattern(model,pattern):
    met_list = []
    for met in model.metabolites:
        if pattern in met.id:
            met_list.append(met.id)
    return met_list

def flux_based_reactions(model,met_id,only_types=(),ignore_types = (),threshold = 0.,flux_dict=0,growth_symbol='mu'):
    if not flux_dict:
        flux_dict = model.solution.x_dict
    reactions = get_reactions_of_met(model,met_id,only_types=only_types,ignore_types=ignore_types,verbose=False)
    if len(reactions) == 0:
        print('No reactions found for {}'.format(met_id))
        return

    result_dict = {}
    for rxn in reactions:
        result_dict[rxn.id] = {}
        for rxn_met,stoich in rxn.metabolites.items():
            if rxn_met.id == met_id:
                if hasattr(stoich, 'subs'):
                    coeff = float(stoich.subs(growth_symbol,flux_dict['biomass_dilution']))
                else:
                    coeff = stoich
                result_dict[rxn.id]['lb'] = rxn.lower_bound
                result_dict[rxn.id]['ub'] = rxn.upper_bound
                result_dict[rxn.id]['rxn_flux'] = flux_dict[rxn.id]
                result_dict[rxn.id]['met_flux'] = flux_dict[rxn.id]*coeff
                try: result_dict[rxn.id]['reaction'] = rxn.reaction
                except: result_dict[rxn.id]['reaction'] = 'no_reaction'
                break
    df = pd.DataFrame.from_dict(result_dict).T
    return df.loc[df['met_flux'].abs().sort_values(ascending=False).index]

def generate_gene_field(me):
    import cobra
    current_gene_ids = [gene.id for gene in me.genes]
    for met in me.metabolites:
        met_id = met.id
        if isinstance(met, cobrame.TranslatedGene):
            gene_id = met_id.split('_')[1]
            if gene_id and gene_id not in current_gene_ids:
                try:
                    gene = cobra.Gene(gene_id)
                    me.genes.append(gene)
                    print(gene_id)
                except:
                    pass

def solution_summary(me):
    reactions = [rxn.id for rxn in me.reactions]
    summary_df = pd.DataFrame(columns=['lb','ub','flux','formula'],index=reactions)

    for rxn_id in tqdm(reactions):
        rxn = me.reactions.get_by_id(rxn_id)
        summary_df.loc[rxn_id]['lb'] = rxn.lower_bound
        summary_df.loc[rxn_id]['ub'] = rxn.upper_bound
        summary_df.loc[rxn_id]['flux'] = me.solution.x_dict[rxn_id]
        summary_df.loc[rxn_id]['formula'] = rxn.reaction

    return summary_df

def get_flux_for_escher(model,type='m'):
    if type == 'm':
        flux_dict = model.solution.x_dict
    elif type == 'me':
        flux_dict = model.get_metabolic_flux(solution=me.solution)

    return pd.DataFrame.from_dict({'flux':flux_dict})

def get_compartments_of_reaction(r):
    return r.get_compartments()

def get_all_transport_of_model(model):
    transport_reactions = []
    for r in tqdm(model.reactions):
        comps = get_compartments_of_reaction(r)
        if len(comps) > 1:
            transport_reactions.append(r.id)
    return list(set(transport_reactions))


def get_transport_reactions(model,met_id,comps=['e','c']):
    from_met = re.sub('_[a-z]$','_'+comps[0],met_id)
    to_met = re.sub('_[a-z]$','_'+comps[1],met_id)

    if isinstance(model,MEModel):
        reaction_type = ['MetabolicReaction']
    else:
        reaction_type = 0
    prod_rxns = [rxn.id for rxn in get_reactions_of_met(model,to_met,s=1,verbose=0,only_types=reaction_type)]
    cons_rxns = [rxn.id for rxn in get_reactions_of_met(model,from_met,s=-1,verbose=0,only_types=reaction_type)]
    transport_rxn_ids = list(set(prod_rxns)&set(cons_rxns))

    transport_rxns = [model.reactions.get_by_id(rxn_id) for rxn_id in transport_rxn_ids]


    return transport_rxns

def global_mass_balance(model):
    exchange_df = exchange_single_model(model)
    balance = 0
    mass_dict = {}
    for r_id in exchange_df.index:
        r = model.reactions.get_by_id(r_id)
        flux = model.solution.x_dict[r_id] # mmol/gDW h
        mass = 0
        for m in r.metabolites:
            coeff = r.metabolites[m]/1000
            weight = m.formula_weight # g/mmol
            mass += coeff*flux*weight
            balance += mass
        mass_dict[r_id] = mass
    return balance, pd.DataFrame.from_dict({'mass':mass_dict})


def get_gene_annotation(me,model):
    gene_annotation = {}
    for m in tqdm(me.metabolites):
        if isinstance(m,cobrame.TranslatedGene):
            gene_id = m.id.split('_')[-1]
            if hasattr(model.genes,gene_id):
                gene = model.genes.get_by_id(gene_id)
                for r in gene.reactions:
                    gene_annotation[gene_id] = r.subsystem
            else:
                rxns = get_reactions_of_met(me,m.id,verbose=False)
                for r in rxns:
                    if 'formation' in r.id:
                        active_complex = [i.id for i in r.products][0]
                        final_rxns = get_reactions_of_met(me,active_complex,verbose=False)
                        subsystem = list(set([i.id.split('_')[0] for i in final_rxns if active_complex not in i.id]))
                        if subsystem:
                            gene_annotation[gene_id] = subsystem[0]
                        break
    return gene_annotation

def get_final_reactions_of_gene(me,gene_id):
    rxns = get_reactions_of_met(me,'protein_'+gene_id,verbose=False)
    final_rxns = []
    if hasattr(me.reactions,'translocation_'+gene_id):
        translocated_complex = [i.id for i in me.reactions.get_by_id('translocation_'+gene_id).products if 'Membrane' in i.id][0]
        formation_rxn = get_reactions_of_met(me,translocated_complex,verbose=False,s=-1)[0]
    else:
        for r in rxns:
            if 'formation' in r.id:
                formation_rxn = r
                break
    active_complex = [i.id for i in formation_rxn.products if 'biomass' not in i.id][0]
    final_rxns = get_reactions_of_met(me,active_complex,verbose=False)
    return list(set(final_rxns) - set([formation_rxn]))

def get_met_production(model,met_list,flux_responses,x_var,only_types = [],plot=True):

    if plot:
        fig,ax = plt.subplots(int(np.ceil(np.sqrt(len(met_list)))),int(np.floor(np.sqrt(len(met_list)))),figsize=(13,3*int(np.ceil(np.sqrt(len(met_list))))))

        if len(met_list) > 1: 
            ax = ax.flatten()

    if not isinstance(flux_responses,list):
        flux_responses = [flux_responses]

    for idx,met_id in enumerate(tqdm(met_list,position=0,leave=True)):
        met = model.metabolites.get_by_id(met_id)

        if isinstance(met,cobrame.Metabolite) and not only_types:
            only_types = ['MetabolicReaction','MEReaction']

        for flux_df in flux_responses:
            met_rate = []
            uptake_rate = []
            for case in flux_df.columns:
                df = flux_based_reactions(model,met_id,flux_dict=flux_df[case].to_dict(),only_types=only_types)
                met_rate.append(df[df.met_flux>0]['met_flux'].sum())
                uptake_rate.append(flux_df.abs()[case][x_var])
            yield met_id,met_rate,uptake_rate

            if plot:
                ax_i = ax[idx] if len(met_list) > 1  else ax
                ax_i.plot(uptake_rate,met_rate,'-o')
                ax_i.set_xlabel(x_var)
                ax_i.set_ylabel('mmol {}/gDW/h'.format(met_id))
                ax_i.set_title(met_id)
        if plot: fig.tight_layout()


def get_compartment_transport(model,comps,MEmodel=False):
    reactions = []
    comps = set(comps)
    for r in model.reactions:
        if MEmodel and not isinstance(r,cobrame.MetabolicReaction):
            continue
        r_comps = set()
        for m in r.metabolites:
            r_comps.add(m.id[-1])
        if len(r_comps & comps)>0 and len(r_comps) > 1:
            reactions.append(r)
    return reactions

def get_biomass_fractions(model,bof):
    fractions = {}
    fractions['protein'] = 0
    fractions['rna'] = 0
    fractions['dna'] = 0
    fractions['lipid'] = 0
    fractions['other'] = 0
    
    
    aminoacid_exp = '^[a-z]{3}__L_c$'
    dna_exp = '^d[a,t,c,g]tp_c$'
    rna_exp = '^[a,t,c,g,u][t,m]p_c$'
    lipid_exp = '^.*_BS_c$'
    for m,coeff in bof.metabolites.items():
        if coeff > 0 or 'h2o' in m.id: continue
        if re.search(aminoacid_exp,m.id):
            fractions['protein'] -= coeff*m.formula_weight
        elif re.search(dna_exp,m.id):
            fractions['dna'] -= coeff*m.formula_weight
        elif re.search(rna_exp,m.id):
            if m.id == 'atp_c':
                coeff = coeff + bof.metabolites[model.metabolites.adp_c]
            fractions['rna'] -= coeff*m.formula_weight
        elif re.search(lipid_exp,m.id):
            fractions['lipid'] -= coeff*m.formula_weight
        else:
            fractions['other'] -= coeff*m.formula_weight        
    return pd.DataFrame.from_dict({'fraction':fractions})
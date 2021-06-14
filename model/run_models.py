import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import prep_for_model_runs as prep
import build_models as build
import modify_contact as mod
import model_params_class as mp
import calc_summary_stat as summ



"""
Run all of the models of interest

params
------
base_dir: String
    The path for the directory to which results should be saved

params: object of class ModelParams
    the model parameters for this model run

days: Int
    the number of days to run the model for

group_size_data: DataFrame (???)
    the DataFrame that contains the size of each population sub-group

contact_data_pre_SIP: k*k contact matrix for pre-SIP order

contact_data_post_SIP: k*k contact matrix for post-SIP order

prison_peak_date: Int
    the date on which the prison infection rate will peak

returns
-------
summary_stats_original: a DataFrame of summary statistics

infection_rates: a DataFrame of all the infection rates for each version of the
model 

pop_size_df: the population sizes used for inputs to each of the models 
"""

def run_model_simple(output_dir, params, days,
               group_size_data,
               contact_data_pre_SIP,
               contact_data_post_SIP):
    output_name = params.get_name()
    output_path =os.path.join(output_dir, output_name)
     
    #get sizes and initial values     
    group_size_data = prep.process_group_size(group_size_data, params.initial_infection_multiplier)        
    pop_size = prep.combine_groups(group_size_data)
      
    #conctact_matrix_list 
    cm_pre = []
    cm_post = []

    pre_lockdown, post_lockdown, group_size_ls = mod.create_matrix(contact_data_pre_SIP,
                                                               contact_data_post_SIP,
                                                               group_size_data)
    infection_rates = []
    s_rates = []
    pop_sizes = {}
    #run the original model
    S_df, I_df, _ = build.build_model(
        group_size_data, days, 
        contact_data_pre_SIP.values,
        contact_data_post_SIP.values, params)
    
    pop_sizes[output_name] = pop_size
    
    cm_pre.append(prep.add_contact_matrix(contact_data_pre_SIP, output_name, 'original'))
    cm_post.append(prep.add_contact_matrix(contact_data_post_SIP, output_name, 'original'))            
                    
    summary_stats_original = summ.calculate_peak_infections(S_df, I_df, pop_size, "original", output_path)
    peak_days_original = summary_stats_original.loc['Total_w_Police', 'days_peak']
    
    for d in [I_df, S_df]:
        d['model_name'] = 'original'
        d['name'] = output_name
    infection_rates.append(I_df)
    s_rates.append(S_df)
       
    #Turn map of population sizes into data frame 
    pop_series = [values.rename(k) for k, values in pop_sizes.items()]
    pop_size_df = pd.concat(pop_series, axis = 1).transpose()
    summary_stats_original['model_tag'] = output_name
    return summary_stats_original, pd.concat(infection_rates), pd.concat(s_rates), pop_size_df



"""
Run all of the models of interest

params
------
base_dir: String
    The path for the directory to which results should be saved

params: object of class ModelParams
    the model parameters for this model run

days: Int
    the number of days to run the model for

group_size_data: DataFrame (???)
    the DataFrame that contains the size of each population sub-group

contact_data_pre_SIP: k*k contact matrix for pre-SIP order

contact_data_post_SIP: k*k contact matrix for post-SIP order

prison_peak_date: Int
    the date on which the prison infection rate will peak

returns
-------
summary_stats_original: a DataFrame of summary statistics

infection_rates: a DataFrame of all the infection rates for each version of the
model 

pop_size_df: the population sizes used for inputs to each of the models 
"""

def run_models(output_dir, params, days,
               group_size_data,
               contact_data_pre_SIP,
               contact_data_post_SIP):
    output_name = params.get_name()
    output_path =os.path.join(output_dir, output_name)
     
    #get sizes and initial values     
    group_size_data = prep.process_group_size(group_size_data, params.initial_infection_multiplier)        
    pop_size = prep.combine_groups(group_size_data)
      
    #conctact_matrix_list 
    cm_pre = []
    cm_post = []

    pre_lockdown, post_lockdown, group_size_ls = mod.create_matrix(contact_data_pre_SIP,
                                                               contact_data_post_SIP,
                                                               group_size_data)
    infection_rates = []
    s_rates = []
    pop_sizes = {}
    #run the original model
    S_df, I_df, _ = build.build_model(
        group_size_data, days, 
        contact_data_pre_SIP.values,
        contact_data_post_SIP.values, params)
    
    pop_sizes['no_policy_original'] = pop_size
    
    cm_pre.append(prep.add_contact_matrix(contact_data_pre_SIP, output_name, 'original'))
    cm_post.append(prep.add_contact_matrix(contact_data_post_SIP, output_name, 'original'))            
                    
    summary_stats_original = summ.calculate_peak_infections(S_df, I_df, pop_size, "original", output_path)
    peak_days_original = summary_stats_original.loc['Total_w_Police', 'days_peak']
    
    for d in [I_df, S_df]:
        d['model_name'] = 'original'
        d['Policy_Lever'] = 'no_policy'
        d['name'] = output_name
    infection_rates.append(I_df)
    s_rates.append(S_df)
    for model_name in pre_lockdown.keys():
        
        pre_contact = pre_lockdown[model_name]
        post_contact = post_lockdown[model_name]
        size_df = group_size_ls[model_name]
                   
        cm_pre.append(prep.add_contact_matrix(contact_data_pre_SIP, output_name, model_name))
        cm_post.append(prep.add_contact_matrix(contact_data_post_SIP, output_name, model_name))  
        
        S_df2, I_df2, _ = build.build_model(
            size_df, days,
            pre_contact.values, post_contact.values, params)
        pop_size = prep.combine_groups(size_df)
        summary_stats = summ.calculate_peak_infections(S_df2, I_df2, pop_size, model_name,
                                                  output_path, peak_days_original)

        for d in [I_df2, S_df2]:
            d['model_name'] = model_name
            d['Policy_Lever'] = 'no_policy'
            d['name'] = output_name
        infection_rates.append(I_df2)
        s_rates.append(S_df2)
        pop_sizes[model_name] = pop_size
        summary_stats_original = summary_stats_original.append(summary_stats)
        
    #Turn map of population sizes into data frame 
    pop_series = [values.rename(k) for k, values in pop_sizes.items()]
    pop_size_df = pd.concat(pop_series, axis = 1).transpose()
    pop_size_df['name'] = output_name
    summary_stats_original['model_tag'] = output_name
        
    return summary_stats_original, pd.concat(infection_rates), pd.concat(s_rates), pop_size_df



# Run the policy interventions. Has if statement for both kinds of policies.
# If the optional params are included is policy lever 1. 
def run_policy_intervention(policy_name,
    output_dir, params, days,
    group_size_data,
    contact_data_pre_SIP, contact_data_post_SIP,
    jail_release_shrink, jail_release_date, policy2_params = None):
    
    output_name = policy_name
    output_path =os.path.join(output_dir, output_name)
    
    #get sizes and initial values     
    group_size_data = prep.process_group_size(group_size_data, params.initial_infection_multiplier)        
    pop_size = prep.combine_groups(group_size_data)
    
    pop_sizes = {}
    
    if policy2_params == None:
        #run the original model
        print("policy_lever_1")
        S_df, I_df, _ = build.build_model_p1(
            group_size_data, days, params.sip_start_date,
            contact_data_pre_SIP.values, contact_data_post_SIP.values,
            params.transmission_rate,params.post_sip_transmission_rate, params.prison_infection_rate,
            params.prison_peak_date, jail_release_shrink, jail_release_date)

# def build_model_p1(group_size_data, TIME, SIP_DATE, contact_matrix1, contact_matrix2,
#                 transmission_rate, post_sip_transmission_rate, prison_peak_rate, prison_peak_date, jail_release_shrink,jail_release_date):
    else:
        print("policy_lever_2")
        
# def build_model_p2(group_size_data, TIME, SIP_DATE, contact_matrix1, contact_matrix2,
#                 transmission_rate, post_sip_transmission_rate, prison_peak_rate, prison_peak_date,
#                 jail_release_date, policy2_params):
        S_df, I_df, _ = build.build_model_p2(
            group_size_data, days,
            contact_data_pre_SIP.values, contact_data_post_SIP.values, params, policy2_params)
        
    pop_sizes[f'{policy_name}_original'] = pop_size
         
    summary_stats_original = summ.calculate_peak_infections(S_df, I_df, pop_size, 'original', output_path)
    peak_days_original = summary_stats_original.loc['Total_w_Police', 'days_peak']
    for d in [I_df, S_df]:
        d['model_name'] = 'original'
        d['Policy_Lever'] = policy_name
        d['name'] = output_name

    #Turn map of population sizes into data frame 
    pop_series = [values.rename(k) for k, values in pop_sizes.items()]
    pop_size_df = pd.concat(pop_series, axis = 1).transpose()
    pop_size_df['name'] = params.get_name()
    summary_stats_original['model_tag'] = output_name
        
    return summary_stats_original, I_df, S_df, pop_size_df
    
    

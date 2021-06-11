# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import prep_for_model_runs as prep
import build_models as build
import modify_contact as mod
import model_params_class as mp
import calc_summary_stat as summ


"""Build model for policy intervention 1
This function is the same as build_model in build_models.py with the exception 
of the additional parameter, group_size_p1. This is the group_df_p1 data.frame
from generate_matrix.py that reduces the prison release population a certain
number of days after the SIP_DATE. 
"""
def build_model_p1(group_size_data, TIME, SIP_DATE, contact_matrix1, contact_matrix2,
                transmission_rate, prison_peak_rate, prison_peak_date, group_size_p1,
                jail_release_date):
    
    Group_Names, Group_Size, initial_sizes, recovery_rates = build.build_initial_values(group_size_data)
    
    susceptible_rows = []
    infected_rows = []
    recovered_rows = []
    lambda_matrix = contact_matrix1 * transmission_rate
 
    S_t, I_t, R_t = initial_sizes
    susceptible_rows.append(S_t)
    infected_rows.append(I_t)
    recovered_rows.append(R_t)
    
    white_prison_i = np.where(Group_Names == 'White_Prison')
    black_prison_i = np.where(Group_Names == 'Black_Prison')
    
    k1, k2 = prison_rate_build(
        Group_Size, prison_peak_date, white_prison_i, black_prison_i, prison_peak_rate)
  
    #represents new infections per day
    delta_i = [R_t]
    for i in range(0, TIME):
        
        if i == SIP_DATE - 1:
            lambda_matrix = contact_matrix2 * transmission_rate
            
        if i == jail_release_date - 1:
            Group_Size = group_size_pol1['Group_Size'].values
            # NOT SURE IF THIS IS CORRECT COLUMN BECAUSE I'M NOT SURE IF WE'RE
            # USING THE OLD FORMAT OF THE GROUP SIZE DATA OR NOT
            k1, k2 = prison_rate_build(Group_Size, prison_peak_date, white_prison_i,
                                       black_prison_i, prison_peak_rate)
            
        # multiplying k*k contact matrix * k*1 vetor of proportion of group infected
        #l is a vector with length k 
        l = np.squeeze(np.asarray(np.dot(lambda_matrix, I_t/Group_Size)))
        
        #this is the number of new infections 
        contacts = l * S_t #force of infection * number Susceptible by group
        delta_i.append(contacts)
        
        I_14 = R_t[0]
        if i >= 14:
            I_14 = delta_i[i-14]
        
        dSdt = - contacts 
        dIdt = contacts - recovery_rates * I_14       
        dRdt = recovery_rates * I_14

        S_t = S_t + dSdt   
        I_t = I_t + dIdt
        R_t = R_t + dRdt
        
        if i <= prison_peak_date:
            I_t[white_prison_i] = np.exp(i*k1)
            I_t[black_prison_i] = np.exp(i*k2)
            S_t[white_prison_i] = Group_Size[white_prison_i] - np.exp(i*k1)
            dSdt[black_prison_i] = Group_Size[black_prison_i] - np.exp(i*k2)
            # Should this be S_t?
        
        susceptible_rows.append(S_t)
        infected_rows.append(I_t)
        recovered_rows.append(R_t)
    s = pd.DataFrame(susceptible_rows, columns=Group_Names)
    i = pd.DataFrame(infected_rows, columns=Group_Names)
    r = pd.DataFrame(recovered_rows, columns=Group_Names)
    return s,i,r

"""
Run all of the models of interest for policy lever 1

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

eq_group_size_data: DataFrame
    the DataFrame where some population sub-groups are equal (???)

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

def run_models_p1(base_dir, params, days,
               group_size_data,
               eq_group_size_data,
               contact_data_pre_SIP,
               contact_data_post_SIP,
               prison_peak_date, group_size_p1, 
               jail_release_date):
    output_name = params.get_name()
    output_path =os.path.join(base_dir, output_name)
     
    #get sizes and initial values     
    group_size_data = prep.process_group_size(group_size_data, params.initial_infection_multiplier)        
    pop_size = prep.combine_groups(group_size_data)
    #NOTE: NOT SURE IF WE SHOULD ALSO RUN THESE ON GROUP_SIZE_P1
   
    #equal forced labour  
    eq_group_size_data = prep.process_group_size(eq_group_size_data, params.initial_infection_multiplier) 
    
    #conctact_matrix_list 
    cm_pre = []
    cm_post = []
    pre_lockdown, post_lockdown, group_size_ls = mod.create_matrix(contact_data_pre_SIP,
                                                               contact_data_post_SIP,
                                                               group_size_data,
                                                               eq_group_size_data)
    #NOTE: DO WE NEED TO MODIFY ABOVE CODE FOR GROUP_SIZE_P1? 
    infection_rates = []
    pop_sizes = {}
    
    #run the original model
    S_df, I_df, _ = build_model_p1(group_size_data, days,
                                params.sip_start_date,
                                contact_data_pre_SIP.values,
                                contact_data_post_SIP.values,
                                params.transmission_rate,
                                params.prison_infection_rate, 
                                prison_peak_date, group_size_p1,
                                jail_release_date)
    pop_sizes['original'] = pop_size
    
    cm_pre.append(prep.add_contact_matrix(contact_data_pre_SIP, output_name, 'original'))
    cm_post.append(prep.add_contact_matrix(contact_data_post_SIP, output_name, 'original'))            
                    
    summary_stats_original = summ.calculate_peak_infections(S_df, I_df, pop_size, "original", output_path)
    peak_days_original = summary_stats_original.loc['Total_w_Police', 'days_peak']
    I_df['model_name'] = 'original'
    I_df['policy'] = 'policy1'
    I_df['name'] = output_name
    infection_rates.append(I_df)
   
    #Turn map of population sizes into data frame 
    pop_series = [values.rename(k) for k, values in pop_sizes.items()]
    pop_size_df = pd.concat(pop_series, axis = 1).transpose()
    pop_size_df['name'] = output_name
    summary_stats_original['model_tag'] = output_name
        
    return summary_stats_original, pd.concat(infection_rates), pop_size_df
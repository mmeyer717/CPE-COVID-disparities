
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import prep_for_model_runs as prep
import build_models as build
import modify_contact as mod
import model_params_class as mp
import calc_summary_stat as summ


"""Build model for policy intervention 2
This function is the same as build_model in build_models.py with the exception 
of the additional parameter that it requires an updated infection rate for prison churn. 
"""
def build_model_p2(group_size_data, TIME, contact_matrix1, contact_matrix2, params, policy2_params):
    

    
    Group_Names, Group_Size, initial_sizes, recovery_rates = build_initial_values(group_size_data)
    
    susceptible_rows = []
    infected_rows = []
    recovered_rows = []
    lambda_matrix = contact_matrix1 * params.transmission_rate
 
    S_t, I_t, R_t = initial_sizes
    susceptible_rows.append(S_t)
    infected_rows.append(I_t)
    recovered_rows.append(R_t)
    
    white_prison_i = np.where(Group_Names == 'White_Prison')
    black_prison_i = np.where(Group_Names == 'Black_Prison')
    
    k1, k2 = prison_rate_build(
        Group_Size, params.prison_peak_date, white_prison_i, black_prison_i, params.prison_infection_rate)
  
    #represents new infections per day
    delta_i = [R_t]
    for i in range(0, TIME):
        
        if i == params.sip_start_date - 1:
            lambda_matrix = contact_matrix2 * params.post_sip_transmission_rate
             
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
        
        if i <= params.prison_peak_date:
            if i <= params.sip_start_date - 1:
                I_t[white_prison_i] = np.exp(i*k1)
                I_t[black_prison_i] = np.exp(i*k2)
                S_t[white_prison_i] = Group_Size[white_prison_i] - np.exp(i*k1)
                S_t[black_prison_i] = Group_Size[black_prison_i] - np.exp(i*k2)
            else:
                print(f'{i}: now reducing prison rates')
                I_t[white_prison_i] = (policy2_params.prison_sip_i_white + policy2_params.jail_sip_i_white)*Group_Size[white_prison_i]
                I_t[black_prison_i] = (policy2_params.prison_sip_i_black + policy2_params.jail_sip_i_black)*Group_Size[black_prison_i]
                S_t[white_prison_i] = Group_Size[white_prison_i] - (policy2_params.prison_sip_i_white + policy2_params.jail_sip_i_white)*Group_Size[white_prison_i]
                S_t[black_prison_i] = Group_Size[black_prison_i] - (policy2_params.prison_sip_i_black + policy2_params.jail_sip_i_black)*Group_Size[black_prison_i]
        
        susceptible_rows.append(S_t)
        infected_rows.append(I_t)
        recovered_rows.append(R_t)
    s = pd.DataFrame(susceptible_rows, columns=Group_Names)

    i = pd.DataFrame(infected_rows, columns=Group_Names)
    r = pd.DataFrame(recovered_rows, columns=Group_Names)
    
    s['Day'] = s.index
    i['Day'] = i.index
    r['Day'] = r.index
    return s,i,r
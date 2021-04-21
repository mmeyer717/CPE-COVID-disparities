import numpy as np
import pandas as pd

"""
Parse groups size DataFrame into the initial values required for the model

Params
------
groupsize_df: DataFrame 
    Includes columns for 'Population_Size', 'Recovery_Rate',
    'Initial_Infection_Rate' for each of the groups.

return
-------
gn: array with the group names (index of the input dataframe)

gs: array with the group sizes 

y0: tuple with Susceptible, Infection, Recovered at time zero. The length of each
    array is equal to the number of groups
    Initially no one has recovered an all but the initially infected are susceptible
"""
def build_initial_values(groupsize_df):

    gn = groupsize_df.index.values

    groups = len(gn) #number of groups 
    
    # proportion of group i that contacts group b
    gs = groupsize_df['Population_Size'].values

    recovery_rates = groupsize_df['Recovery_Rate'].values 

    # Initial number of infected and recovered individuals for each group 
    I0, R0 = groupsize_df['Initial_Infection_Rate'].values*gs, np.zeros(groups)
    # Everyone else, S0, is susceptible to infection initially.
    S0 = gs - I0

    y0 = S0, I0, R0
    
    return gn, gs, y0, recovery_rates


"""
Calculate the k value for exponential growth in prison rate infections for the 
two prison population sizes. 

params
------
group_size: pandas series with sizes of each group

prison_peak_date: Int 
    number of days until prison infection rate hits its peak

prison_index1: Int
    index where the prison population for whites is in group size vector

prison_index2: Int
    index where prison population for blacks is in group size vector

prison_peak_rate: Float
    the rate that the prison infection rate will end at by the prison_peak_date

return
------
k1: The exponential constant of increase for the first group (Whites)
    I_1 = I^(k1 * i)  where i is in the days since the model started.

k2: The exponential constant of increase for the second group  (Blacks)

"""
def prison_rate_build(
    group_size, prison_peak_date, prison_index1, prison_index2, prison_peak_rate):
    t1 = group_size[prison_index1]*prison_peak_rate
    t2 = group_size[prison_index2]*prison_peak_rate
    k_1 = np.log((t1))/prison_peak_date
    k_2 = np.log((t2))/prison_peak_date
    return k_1, k_2



"""Build model for policy intervention 1
This function is the same as build_model in build_models.py with the exception 
of the additional parameter, group_size_p1. This is the group_df_p1 data.frame
from generate_matrix.py that reduces the prison release population a certain
number of days after the SIP_DATE. 
"""
def build_model_p1(group_size_data, TIME, SIP_DATE, contact_matrix1, contact_matrix2,
                transmission_rate, post_sip_transmission_rate, prison_peak_rate, prison_peak_date, jail_release_shrink,jail_release_date):
    
    Group_Names, Group_Size, initial_sizes, recovery_rates = build_initial_values(group_size_data)
    
    susceptible_rows = []
    infected_rows = []
    recovered_rows = []
    prison_infection = {"Black":[], "White":[]}
    lambda_matrix = contact_matrix1 * transmission_rate
 
    S_t, I_t, R_t = initial_sizes
    susceptible_rows.append(S_t)
    infected_rows.append(I_t)
    recovered_rows.append(R_t)
    
    white_prison_i = np.where(Group_Names == 'White_Prison')
    black_prison_i = np.where(Group_Names == 'Black_Prison')
    
    k1, k2 = prison_rate_build(
        Group_Size, prison_peak_date, white_prison_i, black_prison_i, prison_peak_rate)
    # for use in shrinking jail size 
    jail_release_shrink_by_day = jail_release_shrink/(jail_release_date - SIP_DATE)
    orig_prison_pop_white = Group_Size[white_prison_i]
    orig_prison_pop_black = Group_Size[black_prison_i]
    JAIL_OF_CORRECTIONS = 27296/(27296+1704) # fraction of jail/prison releases that are jail releases;
# comes from Calculations tab, cells C24 and C27  
    #represents new infections per day
    delta_i = [R_t]
    days_since_lockdown = 1
    for i in range(0, TIME):
        
        if i == SIP_DATE - 1:
            lambda_matrix = contact_matrix2 * post_sip_transmission_rate
        
        if i >= SIP_DATE and i <= jail_release_date - 1: #if date after SIP and before/on final jail release date
            Group_Size[white_prison_i] = orig_prison_pop_white*(1-(
                JAIL_OF_CORRECTIONS*jail_release_shrink_by_day*days_since_lockdown))
            Group_Size[black_prison_i] = orig_prison_pop_black*(1-(
                JAIL_OF_CORRECTIONS*jail_release_shrink_by_day*days_since_lockdown))
            #don't recalculate k1, k2
            k1, k2 = prison_rate_build(Group_Size, prison_peak_date, white_prison_i,
                                       black_prison_i, prison_peak_rate)
            days_since_lockdown += 1
            
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
        # I just copied and pasted the code from no lever; no modifications yet
        if i <= prison_peak_date:
            prison_infection["White"].append(np.exp(i*k1)) # number of white prisoners newly infected on day i
            prison_infection["Black"].append(np.exp(i*k2)) # number of black prisoners newly infected on day i
            S_t[white_prison_i] = Group_Size[white_prison_i] - np.exp(i*k1)
            S_t[black_prison_i] = Group_Size[black_prison_i] - np.exp(i*k2)
        elif i > prison_peak_date:
            prison_infection["White"].append(prison_infection["White"][-1])
            prison_infection["Black"].append(prison_infection["Black"][-1])
        if i > 9: # make sure only the last ten days are kept in the infections lists
            prison_infection["White"].pop(0) 
            prison_infection["Black"].pop(0)
        I_t[white_prison_i] = sum(prison_infection["White"])
        I_t[black_prison_i] = sum(prison_infection["Black"])   
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

"""
Build a model for new infections each day given initial group size, contacts rates, 
transmision rates etc for k groups. Specifically, output a DataFrame with the 
daily number of infections, recovered, susceptible 

params
--------
group_size_data: array of length k
    population size for each population sub-group with names for each sub-group

TIME: Int
    total number of days to model
    
SIP_DATE: Int 
    number of days until shelter in place order -- will use contact_matrix1 for
    contact rates until this number is reached
    
contact_matrix1: k*k matrix
    expected number of contacts that each group makes with every other group 
    BEFORE the shelter in place order
    
contact_matrix2: k*k matrix
    expected number of contacts that each group makes with every other group 
    AFTER the shelter in place order     

transmission_rate: float
    scalar, represents the likelihood of being infected by each contact
    
prison_peak_rate: Float
    the rate that the prison infection rate will end at by the prison_peak_date

prison_peak_date: Int 
    number of days until prison infection rate hits its peak


return
------
susceptible: k*TIME DataFrame where each column is a sub-group, and each row is
    the number of susceptible people on that day for each group
infected: DataFrame with daily number of active infections for each sub-group and 
    each day 
recovered_rows: DataFrame with number of recovered inidividuals in each sub-group 
    on each day. 

Note: the sum of cells in each corresponding cell in the three DataFrames is the
number of total people in that group. 
"""

def build_model(group_size_data, TIME, contact_matrix1, contact_matrix2,
                params):
    
    Group_Names, Group_Size, initial_sizes, recovery_rates = build_initial_values(group_size_data)
    
    susceptible_rows = []
    infected_rows = []
    recovered_rows = []
    # I use the below dictionary to keep track of daily infections in prison by race
    prison_infection = {"Black":[], "White":[]}
    lambda_matrix = contact_matrix1 * params.transmission_rate
 
    S_t, I_t, R_t = initial_sizes
    susceptible_rows.append(S_t)
    infected_rows.append(I_t)
    recovered_rows.append(R_t)
    
    white_prison_i = np.where(Group_Names == 'White_Prison')
    black_prison_i = np.where(Group_Names == 'Black_Prison')
    
    k1, k2 = prison_rate_build(
        Group_Size, params.prison_peak_date, 
        white_prison_i, black_prison_i,params.prison_infection_rate)
  
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
            #KEEP TRACK OF RELEASES
            prison_infection["White"].append(np.exp(i*k1)) # number of white prisoners newly infected on day i
            prison_infection["Black"].append(np.exp(i*k2)) # number of black prisoners newly infected on day i
            S_t[white_prison_i] = Group_Size[white_prison_i] - np.exp(i*k1)
            S_t[black_prison_i] = Group_Size[black_prison_i] - np.exp(i*k2)
        elif i > params.prison_peak_date:
            prison_infection["White"].append(prison_infection["White"][-1])
            prison_infection["Black"].append(prison_infection["Black"][-1])
        if i > 9: # make sure only the last ten days are kept in the infections lists
            prison_infection["White"].pop(0) 
            prison_infection["Black"].pop(0)
        I_t[white_prison_i] = sum(prison_infection["White"])
        I_t[black_prison_i] = sum(prison_infection["Black"])
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

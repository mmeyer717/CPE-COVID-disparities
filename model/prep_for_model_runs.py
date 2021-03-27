import os

"""
Do initial processing on the group size DataFrame

params
------
group_size_data: series with group size information

initial_multiplier: 

returns
-------
the population size, initial infection rate, and recovery rate for each population subgroup
"""

def process_group_size(group_size_data, initial_multiplier):
    group_size_data['Recovery_Rate'] = 1
    #set so people in prisons don't recover 
    group_size_data.loc['White_Prison', 'Recovery_Rate'] = 0.0 
    #set so people in prisons don't recover 
    group_size_data.loc['Black_Prison', 'Recovery_Rate'] = 0.0 
    group_size_data['Initial_Infection_Rate'] = group_size_data['Initial_Infection_Rate'] * initial_multiplier
    group_size_data.loc['White_Prison', 'Initial_Infection_Rate'] = 0.0 #set prision start rate
    group_size_data.loc['Black_Prison', 'Initial_Infection_Rate'] = 0.0
    return group_size_data[['Population_Size', 'Initial_Infection_Rate', 'Recovery_Rate']]


"""
Combine various population sub-groups' sizes

param
-----
group_size_data: a DataFrame containing the population sizes of each population 
    sub-group

returns
-------
pop_size: group size DataFrame with Total_Residents column and
total resident population size with police (Total_w_Police)
"""
def combine_groups(group_size_data):
    pop_size = group_size_data.copy(True)['Population_Size']
    if 'Black_Forced_Labour' in group_size_data.index:
        pop_size['Total_Residents'] = pop_size['Black'] + \
            pop_size['Black_Forced_Labour'] +  pop_size['White'] + \
                pop_size['White_Forced_Labour']        
    else:
        pop_size['Total_Residents'] = pop_size['Black'] +  pop_size['White']

    pop_size['Total_w_Police'] = pop_size['Total_Residents'] + \
        pop_size['Black_Police'] + pop_size['White_Police']
    return pop_size

"""
Create directories for storing output. 
"""
def create_dirs(base_dir, output_name):
    output_path =os.path.join(base_dir, output_name)
    
     #create output directory
    try:
        os.mkdir(output_path)
    except OSError:
        print ("Creation of the directory %s failed" % output_path)
    else:
        print ("Successfully created the directory %s " % output_path)
        
    #where we output all the contact matrix    
    contact_dir = os.path.join(output_path, 'contacts')
    
    try:
        os.mkdir(contact_dir)
    except OSError:
        print ("Creation of the directory %s failed" % contact_dir)
    else:
        print ("Successfully created the directory %s " % contact_dir)
        
    return output_path, contact_dir


"""These functions are just to parse the input as it is given 
"""
#From the directy get the path for contact matrices without equalizing forced labour groups 
def get_GROUP_SIZE_PATH(pgrp):
    return "GROUP_SIZE_pgrp" + str(pgrp) + "_neq.csv"

#Group size paths with population numbers for equal percent of forced labour
def get_GROUP_SIZE_EQ_PATH(pgrp):
    return "GROUP_SIZE_pgrp" + str(pgrp) + "_eq.csv"

def get_CONTACT_MATRIX_PRE_LOCKDOWN(pc, pgrp):
    return "Contact_Matrix_Pre_SIP_pc%s_pgrp%s.csv"%(pc, pgrp)

def get_CONTACT_MATRIX_POST_LOCKDOWN(pc, pgrp):
    return "Contact_Matrix_Post_SIP_pc%s_pgrp%s.csv"%(pc, pgrp)

"""Read the contact matrix and set index accordingly

params
------
df: DataFrame
    Contact matrix dataframe 

output_name: string
    the parameters going in

model_id: int
    ???

returns
-------
a version of the original DataFrame where the indices are now the model_tag
and model_id (model)
"""
def add_contact_matrix(df, output_name, model_id):
    df2 = df.reset_index()
    df2['model_tag'] = output_name
    df2['model'] = model_id
    df2 = df2.set_index(['model_tag', 'model'])
    return df2

    
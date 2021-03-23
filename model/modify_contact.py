"""Create a copy of the contact matrix where no residents come in contact with 
police

params
-----
c: the k*k contact matrix (k = the number of sub-population groups)

returns
-------
a version of the contact matrix where the contact rate between resident 
sub-population groups and the police is 0
"""
def no_police_contact(c):
    c2 = c.copy(True)
    
    #white police 
    c2.loc['Black', 'White_Police'] = 0
    c2.loc['Black_Forced_Labour', 'White_Police'] = 0
    c2.loc['White_Forced_Labour', 'White_Police'] = 0.0
    c2.loc['White', 'White_Police'] = 0
    
    #black police
    c2.loc['Black', 'Black_Police'] = 0
    c2.loc['Black_Forced_Labour', 'Black_Police'] = 0
    c2.loc['White_Forced_Labour', 'Black_Police'] = 0.0
    c2.loc['White', 'Black_Police'] = 0
    return c2


"""Create copy of contact matrix where there is no difference in prison churn 
rate

param
-----
c: the k*k contact matrix

returns
-------
a version of the contact matrix where no resident sub-population groups come in 
contact with any of the prison sub-population groups, which implies no difference
in prison churn
"""

def no_prison_churn(c):
    c2 = c.copy(True)
   
    c2.loc['Black', 'White_Prison'] = 0
    c2.loc['Black_Forced_Labour', 'White_Prison'] = 0
    c2.loc['White', 'White_Prison'] = 0.0
    c2.loc['White_Forced_Labour', 'White_Prison'] = 0.0
    
    c2.loc['Black', 'Black_Prison'] = 0
    c2.loc['Black_Forced_Labour', 'Black_Prison'] = 0
    c2.loc['White', 'Black_Prison'] = 0.0
    c2.loc['White_Forced_Labour', 'Black_Prison'] = 0.0
    return c2


"""Create a copy of contact matrix where Black and White prison churn rates are
equal

param
-----
c: the k*k contact matrix

returns
-------
a version of the contact matrix where the Black prison churn rate is set equal 
to the White prison churn rate
"""
def equalize_jail(c):
    c2 = c.copy(True)
    white_prison_churn = c.loc['White', 'White_Prison']
   
    c2.loc['Black', 'Black_Prison'] = white_prison_churn
    c2.loc['Black_Forced_Labour', 'Black_Prison'] = white_prison_churn
    return c2


"""Create a copy of contact matrix where Black/White police contact is equal

param
-----
c: the k*k contact matrix

returns
-------
a version of the contact matrix where residents' contact rates with police of 
the other race are set to be equal across both races, and residents' contact 
rates with police of the same race are set to be equal across both races
"""
def equalize_police(c):
    c2 = c.copy(True)
    white_op_police_contact = c.loc['White', 'Black_Police']
    c2.loc['Black', 'White_Police'] = white_op_police_contact
    c2.loc['Black_Forced_Labour', 'White_Police'] = white_op_police_contact
    
    white_police_contact = c.loc['White', 'White_Police']
    c2.loc['Black', 'Black_Police'] = white_police_contact
    c2.loc['Black_Forced_Labour', 'Black_Police'] = white_police_contact
    return c2


"""Add forced labour population size to the non-forced labour population size

param
-----
group_df: a DataFrame of sub-population groups and population sizes

returns
-------
g_df: a DataFrame where the non-police and non-prison population subgroups are
combined into one column and a single column represents Blacks and Whites.
"""
def no_forced_labour(group_df):
    g_df = group_df.copy(True)
    g_df.loc['White', 'Population_Size'] = group_df.loc['White', 'Population_Size'] + \
        group_df.loc['White_Forced_Labour', 'Population_Size']
    g_df.loc['Black', 'Population_Size'] = group_df.loc['Black', 'Population_Size'] + \
        group_df.loc['Black_Forced_Labour', 'Population_Size']
    
    g_df = g_df.drop(['White_Forced_Labour'])
    g_df = g_df.drop(['Black_Forced_Labour'])
    return g_df


"""
Drop the forced_labour columns and rows from the contact matrix

param
-----
contact_matrix_or: the original k*k contact matrix

returns
-------
contact_matrix: a new version of the contact matrix where the forced labour rows
and columns are dropped
"""
def drop_forced_labour(contact_matrix_or):
    contact_matrix = contact_matrix_or.copy(True)
    contact_matrix = contact_matrix.drop(['White_Forced_Labour', 'Black_Forced_Labour'])
 
    contact_matrix = contact_matrix.drop('White_Forced_Labour', axis = 1)
    contact_matrix = contact_matrix.drop('Black_Forced_Labour', axis = 1)
    return contact_matrix


"""
Create pre- and post-shelter-in-place order contact matrices and equalized 
group_size_matrices for various scenarios

params
------
pre_lockdown: the pre-SIP order contact matrix

post_lockdown: the post-SIP order contact matrix

group_size: the original group_size matrix

eq_group_size: ???

returns
-------
pre_lockdown_matrix: a dictionary keyed by the scenario, with vales as the pre-
lockdown matrix

post_lockdown_matrix: same keys, values are post lockdown matrix 

group_size_matrix: same keys, values are the series with the group sizes 

"""
def create_matrix_eq(pre_lockdown, post_lockdown, group_size, eq_group_size):
    pre_lockdown_matrix = {}
    post_lockdown_matrix = {}
    
    group_size_matrix = {}

    #No police only
    without_police_contact_p_SIP = no_police_contact(pre_lockdown)
    without_police_contact = no_police_contact(post_lockdown)
    
    
    pre_lockdown_matrix['no_police'] = without_police_contact_p_SIP
    post_lockdown_matrix['no_police'] = without_police_contact
   
    #No jail/prison only
    without_prison_churn_p_SIP = no_prison_churn(pre_lockdown)
    without_prison_churn = no_prison_churn(post_lockdown)
    
    
    pre_lockdown_matrix['no_prison'] = without_prison_churn_p_SIP
    post_lockdown_matrix['no_prison'] = without_prison_churn
    
    #No police AND no jail/prison
    without_prison_churn_or_police_p_SIP = no_prison_churn(without_police_contact_p_SIP)
    without_prison_churn_or_police = no_prison_churn(without_police_contact)
    
    
    pre_lockdown_matrix['no_prison_or_police'] = without_prison_churn_or_police_p_SIP
    post_lockdown_matrix['no_prison_or_police'] = without_prison_churn_or_police
    
    #Equalize police only
    eq_police_p_sip = equalize_police(pre_lockdown)
    eq_police = equalize_police(post_lockdown)
    
    pre_lockdown_matrix['eq_police'] = eq_police_p_sip
    post_lockdown_matrix['eq_police'] = eq_police
    
    #Equalize jail/prison only
    eq_prison_p_sip = equalize_jail(pre_lockdown)
    eq_prison = equalize_jail(post_lockdown)
    
    pre_lockdown_matrix['eq_prison'] = eq_prison_p_sip
    post_lockdown_matrix['eq_prison'] = eq_prison
    
    #Equalize police AND jail/prison
    equalized_contacts_pre_SIP = equalize_jail(eq_police_p_sip)
    equalized_contacts = equalize_jail(eq_police)
    
    pre_lockdown_matrix['eq_police_prison'] = equalized_contacts_pre_SIP
    post_lockdown_matrix['eq_police_prison'] = equalized_contacts
    
    for k in pre_lockdown_matrix.keys():
        group_size_matrix[k] = group_size
    
    no_forced_labour_df = no_forced_labour(group_size)
    
    #Equalize essential worker only --> has different group size but same as 
    #original contact matrix
    pre_lockdown_matrix['eq_forced_labour'] = pre_lockdown
    post_lockdown_matrix['eq_forced_labour'] = post_lockdown
    group_size_matrix['eq_forced_labour'] = eq_group_size
    
    #No essential workers only (adding this back in) --> diff group size but 
    #same as original contact matrix
    pre_lockdown_matrix['no_forced_labour'] = drop_forced_labour(pre_lockdown)
    post_lockdown_matrix['no_forced_labour'] = drop_forced_labour(post_lockdown)
    group_size_matrix['no_forced_labour'] = no_forced_labour_df
    
    #No police AND no jail/prison AND no essential workers
    pre_lockdown_matrix['no_police_prison_fl'] = drop_forced_labour(without_police_contact_p_SIP)
    post_lockdown_matrix['no_police_prison_fl'] = drop_forced_labour(without_police_contact)
    group_size_matrix['no_police_prison_fl'] = no_forced_labour_df
    
    #Equalize police AND jail/prison AND essential workers
    pre_lockdown_matrix['eq_police_prison_fl'] = equalized_contacts_pre_SIP
    post_lockdown_matrix['eq_police_prison_fl'] = equalized_contacts
    group_size_matrix['eq_police_prison_fl'] = eq_group_size

    return pre_lockdown_matrix, post_lockdown_matrix, group_size_matrix


"""
Create pre- and post-shelter-in-place order contact matrices and equalized 
group_size_matrices for various scenarios. NOT including the equalized group matrices

params
------
pre_lockdown: the pre-SIP order contact matrix

post_lockdown: the post-SIP order contact matrix

group_size: the original group_size matrix


returns
-------
pre_lockdown_matrix: a dictionary keyed by the scenario, with vales as the pre-
lockdown matrix

post_lockdown_matrix: same keys, values are post lockdown matrix 

group_size_matrix: same keys, values are the series with the group sizes 

"""
def create_matrix(pre_lockdown, post_lockdown, group_size):
    pre_lockdown_matrix = {}
    post_lockdown_matrix = {}
    
    group_size_matrix = {}

    #No police only
    without_police_contact_p_SIP = no_police_contact(pre_lockdown)
    without_police_contact = no_police_contact(post_lockdown)
    
    pre_lockdown_matrix['no_police'] = without_police_contact_p_SIP
    post_lockdown_matrix['no_police'] = without_police_contact
   
    #No jail/prison only
    without_prison_churn_p_SIP = no_prison_churn(pre_lockdown)
    without_prison_churn = no_prison_churn(post_lockdown)
    
    
    pre_lockdown_matrix['no_prison'] = without_prison_churn_p_SIP
    post_lockdown_matrix['no_prison'] = without_prison_churn
    
    #No police AND no jail/prison
    without_prison_churn_or_police_p_SIP = no_prison_churn(without_police_contact_p_SIP)
    without_prison_churn_or_police = no_prison_churn(without_police_contact)
    
    
    pre_lockdown_matrix['no_prison_or_police'] = without_prison_churn_or_police_p_SIP
    post_lockdown_matrix['no_prison_or_police'] = without_prison_churn_or_police
    
   
    for k in pre_lockdown_matrix.keys():
        group_size_matrix[k] = group_size
    
    no_forced_labour_df = no_forced_labour(group_size)
     
    #No essential workers only (adding this back in) --> diff group size but 
    #same as original contact matrix
    pre_lockdown_matrix['no_forced_labour'] = drop_forced_labour(pre_lockdown)
    post_lockdown_matrix['no_forced_labour'] = drop_forced_labour(post_lockdown)
    group_size_matrix['no_forced_labour'] = no_forced_labour_df
    
    #No police AND no jail/prison AND no essential workers
    pre_lockdown_matrix['no_police_prison_fl'] = drop_forced_labour(without_police_contact_p_SIP)
    post_lockdown_matrix['no_police_prison_fl'] = drop_forced_labour(without_police_contact)
    group_size_matrix['no_police_prison_fl'] = no_forced_labour_df

    return pre_lockdown_matrix, post_lockdown_matrix, group_size_matrix
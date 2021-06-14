
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import prep_for_model_runs as prep
import model_params_class as mp
import run_models as run
import sys
import generate_matrix as inputs
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None


"""
Old code to draw the susceptable plots with the policy levers
"""
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
def make_susceptable_plot(s_maps, group_sizes, suptitle, DAY_LIMIT = 120, policies = ['no_policy']):
    n_policies = len(policies)
    fig, ax  = plt.subplots(n_policies+1, 1, figsize = (10, 15))
    no_policy = s_maps['no_policy']
    p = group_sizes.loc['no_policy']
#     no_policy['all_blacks'] = no_policy['Black'] + no_policy['Black_Forced_Labour']
#     no_policy['all_whites'] = no_policy['White'] + no_policy['White_Forced_Labour']
#     p['all_whites'] = p['White'] + p['White_Forced_Labour']
#     p['all_blacks'] = p['Black'] + p['Black_Forced_Labour']
   
    colors = {
        'Black_All': 'b',
        'White_All': 'r',
        'Total_Residents': 'g'
    }
    for i in range(0, n_policies):
        policy_name = policies[i]
        policy_df = s_maps[policy_name].iloc[0: DAY_LIMIT]
        df = policy_df[policy_df['model_name'] == 'original']
        
#         df['all_blacks'] = df['Black'] + df['Black_Forced_Labour']
#         df['all_whites'] = df['White'] + df['White_Forced_Labour']
#         p['all_whites'] = p['White'] + p['White_Forced_Labour']
#         p['all_blacks'] = p['Black'] + p['Black_Forced_Labour']
        for group in ['Black_All', 'White_All']:
            color = colors[group]
            ax[i].plot(df['Day'].values, 1 - (df[group]/p[group]), color, alpha=0.85, lw=3, label = group)
            ax[i].set_title('Policy Lever: ' + policy_name)
            ax[i].legend()
            ax[i].set_ylabel('cumulative infection rate')
        ax[n_policies].plot(df['Day'].values, 1 - (df['Total_Residents']/p['Total_Residents']),
                   list(colors.values())[i], alpha=0.85, lw=3, label = policy_name+'total')
        ax[n_policies].legend()
        ax[n_policies].set_title('All residents')
        
    
    fig.suptitle(suptitle)
    plt.savefig(f'figures/{suptitle}.png')
    
def make_prison_shrink_list(orig_pop, shrink_factor, ndays_shrink, j_of_c):
    pop_list = []
    for j in range(ndays_shrink):
        pop_list_item = orig_pop*(1-(j_of_c*shrink_factor*(j+1)))
        pop_list.append(pop_list_item)
    return pop_list

populations = {
        'Black_Forced_Labour': 'g',
        'Black': 'b',
        'Total_Residents': 'k+',
        'White_Forced_Labour' : 'y',
        'White': 'r'
    }
    
prisons = {
        'White_Prison' : 'y+',
        'White_Police' : 'r+',
        'Black_Prison': 'g',
        'Black_Police': 'b'
    }
    
def plot_all_groups(s, group_sizes, policies = ['no_policy'],
                    DAY_LIMIT = 120):
    n_policies = len(policies)
    fig, ax  = plt.subplots(len(policies),2, figsize = (10, 5 * n_policies))
    p = group_sizes.loc['no_policy']
    maps = [populations, prisons]
    for j in range(0, n_policies):
        policy_name = policies[j]
        policy_df = s[policy_name].iloc[0:DAY_LIMIT]
        df = policy_df[policy_df['model_name'] == 'original']
        
        axj = ax if len(policies) == 1 else ax[j]
        for i in [0,1]:
            for group, color in maps[i].items():
                if group in ["Black_Prison", "White_Prison"] and policy_name == "lever1":
                    jail_release_shrink_by_day = 0.4/7 # 7 = jail_release_date - SIP_DATE
                    orig_prison_pop = p[group]
                    JAIL_OF_CORRECTIONS = 27296/(27296+1704)
                    # create vector of group sizes for each day; will need to repeat values
                    shrink_list = make_prison_shrink_list(orig_prison_pop,
                                                                jail_release_shrink_by_day, 7,
                                                                JAIL_OF_CORRECTIONS)
                    print("shrink list")
                    print(shrink_list)
                    gs = [orig_prison_pop for ii in range(14)] + shrink_list + [shrink_list[-1] for iii in range(120-21)]
                    print("gs before shrink")
                    print(gs[13])
                    print("gs after shrink")
                    print(gs[21:24])
                    print("infections on shrink days")
                    print(df[group][14:21])
                    axj[i].plot(df['Day'].values, 1 - (df[group]/gs), color, alpha=0.85, lw=3, label = group)
                else:
                    axj[i].plot(df['Day'].values, 1 - (df[group]/p[group]), color, alpha=0.85, lw=3, label = group)

            axj[i].legend()
            axj[i].set_ylabel('cumulative infection rate')
            axj[i].set_title(policy_name)
            
    return ax
   

def plot_all_groups_infection(infection_maps, group_sizes, policies = ['no_policy'],
                              DAY_LIMIT = 120, use_rate = True):
    n_policies = len(policies)
    fig, ax  = plt.subplots(len(policies),2, figsize = (10, 5 * n_policies))

    
    p = group_sizes.loc['no_policy']
    maps = [populations, prisons]
    for j in range(0, n_policies):
        policy_name = policies[j]
        policy_df = infection_maps[policy_name].iloc[0: DAY_LIMIT]
        df = policy_df[policy_df['model_name'] == 'original']
        axj = ax if len(policies) == 1 else ax[j]
        for i in [0,1]:
            for group, color in maps[i].items():
                if use_rate:
                    axj[i].plot(df['Day'].values, df[group]/p[group], color, alpha=0.85, lw=3, label = group)
                    axj[i].set_ylabel('Daily Infection Rate')

                else:
                    axj[i].plot(df['Day'].values, df[group], color, alpha=0.85, lw=3, label = group)
                    axj[i].set_ylabel('Daily Infection #')

            axj[i].legend()
            
    return ax





def make_susceptable_plot_model_dif(s_maps, group_sizes, suptitle):
    fig, ax  = plt.subplots(2, 3, figsize = (20, 10))
    no_policy = s_maps['no_policy']
    p = group_sizes.loc['no_policy']
    no_policy['all_blacks'] = no_policy['Black'] + no_policy['Black_Forced_Labour']
    no_policy['all_whites'] = no_policy['White'] + no_policy['White_Forced_Labour']
    p['all_whites'] = p['White'] + p['White_Forced_Labour']
    p['all_blacks'] = p['Black'] + p['Black_Forced_Labour']
    policies = ['no_policy', 'lever1', 'lever2']
    colors = {
        'no_policy': 'k',
        'lever1': 'r',
        'lever2': 'b'
    }
    policy_df = s_maps['no_policy']
    df_original = policy_df[policy_df['model_name'] == 'original']
    df_original['all_blacks'] = df_original['Black'] + df_original['Black_Forced_Labour']
    df_original['all_whites'] = df_original['White'] + df_original['White_Forced_Labour']
    
    
    for policy in ['lever1', 'lever2']:
        policy_name = policy
        policy_df = s_maps[policy_name]
        color = colors[policy]
        df = policy_df[policy_df['model_name'] == 'original']
        p = group_sizes.loc[policy]
        df['all_blacks'] = df['Black'] + df['Black_Forced_Labour']
        df['all_whites'] = df['White'] + df['White_Forced_Labour']
        p['all_whites'] = p['White'] + p['White_Forced_Labour']
        p['all_blacks'] = p['Black'] + p['Black_Forced_Labour']
        groups = ['all_blacks', 'all_whites', 'Total_Residents']
        for i in range(0, 3):
            group = groups[i]
            
            diff = 1 - (df_original[group]/p[group]) - (1 - (df[group]/p[group]))
            ax[0,i].plot(df['Day'].values, diff , color, alpha=0.85, lw=3, label = policy_name)
            ax[0, i].set_title(
                'Difference in rates ' + f'for {group}')
            ax[0, i].legend()
            
            diff_total = (1 - df_original[group]) - (1 - df[group])
            ax[1,i].plot(df['Day'].values, diff_total , color, alpha=0.85, lw=3, label = policy_name)
            ax[1,i].set_title(
                'Difference in cum cases ' + f'for {group}')
            ax[1, i].legend()
    fig.suptitle(suptitle)
    plt.savefig('figures/{suptitle}.png')
    

policy_colors = {
        'no_policy': 'b',
        'lever1': 'r',
        'lever2': 'g'
    }
def make_susceptable_plot_racial_dif(s_maps, group_sizes, suptitle, policies = ['no_policy']):
    n_policies = len(policies)
    fig, ax  = plt.subplots(1,1, figsize = (10, 5 * n_policies))
    no_policy = s_maps['no_policy']
    p = group_sizes.loc['no_policy']
    no_policy['all_blacks'] = no_policy['Black'] + no_policy['Black_Forced_Labour']
    no_policy['all_whites'] = no_policy['White'] + no_policy['White_Forced_Labour']
    p['all_whites'] = p['White'] + p['White_Forced_Labour']
    p['all_blacks'] = p['Black'] + p['Black_Forced_Labour']
   

    for i in range(0, n_policies):
        policy_name = policies[i]
        policy_df = s_maps[policy_name].iloc[0:120]
        df = policy_df[policy_df['model_name'] == 'original']
        p = group_sizes.loc[policy_name]
        df['all_blacks'] = df['Black'] + df['Black_Forced_Labour']
        df['all_whites'] = df['White'] + df['White_Forced_Labour']
        p['all_whites'] = p['White'] + p['White_Forced_Labour']
        p['all_blacks'] = p['Black'] + p['Black_Forced_Labour']

        color = policy_colors[policy_name]
        race_diff = 1 - (df['all_blacks']/p['all_blacks']) - (1 - df['all_whites']/p['all_whites'])
        ax.plot(df['Day'].values, race_diff , color, alpha=0.85, lw=3, label = policy_name)
        ax.set_title(f'Black - White Cumulative Cases (Policy Lever: {policy_name})')
        ax.set_ylabel('Difference in cumulative infection rate by race')
        ax.legend()
        
    ax.set_title(suptitle)
    plt.savefig(f'figures/{suptitle}.png')
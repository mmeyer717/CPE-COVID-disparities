import pandas as pd
import numpy as np

"""
Given a DataFrame with the number of susceptible people each day, and the number
of daily infected, compute some summary statistics and write out the
number of active infections to a file. 

Params
------
S: DataFrame
    the number of susceptible people each day in the sub-population

I: DataFrame 
    the number of active infections each day in each sub-population

pop_size: DataFrame
    contains the size of each sub-population group

name: String
    Name of the model, e.g. 'No_Police'

base_dir: String
    Directory to write the results to

days_peak_original = None: 
    The day that the number of infection peaked in the "original" model.

returns
-------
summary statistics DataFrame
"""

def calculate_peak_infections(S, I, pop_size, name, base_dir, days_peak_original = None):
    
    if 'Black_Forced_Labour' not in I.columns:
        I['Black_Forced_Labour'] = 0.0
        S['Black_Forced_Labour'] = 0.0
        pop_size['Black_Forced_Labour'] = 0.0
        
        I['White_Forced_Labour'] = 0.0
        S['White_Forced_Labour'] = 0.0
        pop_size['White_Forced_Labour'] = 0.0
    
    
    I["Total_Residents"] = I['Black'] + I['Black_Forced_Labour'] + I['White'] + I['White_Forced_Labour']
    I['Total_w_Police'] = I["Total_Residents"] + I["Black_Police"] + I["White_Police"]
        
    #Add new columns to susceptible
    S["Total_Residents"] = S['Black'] + S['Black_Forced_Labour'] + S['White'] + S['White_Forced_Labour']
    S['Total_w_Police'] = S["Total_Residents"] + S["Black_Police"] + S["White_Police"]
    
    columns = ['Black', 'Black_Forced_Labour','White', 'White_Forced_Labour',
        'Total_Residents', 'Total_w_Police']
    
    new_sizes = pop_size[['Black', 'Black_Forced_Labour','White', 'White_Forced_Labour',
        'Total_Residents', 'Total_w_Police']]
       
    cumulative_df = pd.DataFrame()
    for c in columns:
        cumulative_df[c] = np.round(pop_size[c] - S[c])
    
    peak_infected = I[columns].max()
    days_peak = I[columns].idxmax()
    if days_peak_original is None:
        days_peak_original = days_peak['Total_w_Police']
        
    before_peak_day = days_peak_original - 2
    if days_peak_original + 2 <= len(cumulative_df):
        after_peak_day = days_peak_original + 2
    else:
        after_peak_day = days_peak_original
    active_total_40_days = I[columns].iloc[40]
    
    cum_total_40_days = cumulative_df[columns].iloc[40]
    cum_rate_40_days = cum_total_40_days/new_sizes
    
    cum_total_100_days = cumulative_df[columns].iloc[120]
    
    #before peak
    cum_before = cumulative_df[columns].iloc[before_peak_day]
    cum_rate_before = cum_before/new_sizes
    active_rate_before = I[columns].iloc[before_peak_day]/new_sizes
    
    #at peak 
    cum_at_peak = np.array([cumulative_df.loc[days_peak[c], c] for c in columns])
    cum_at_peak_rate = cum_at_peak/new_sizes
        
    #after peak
    print("%f" % after_peak_day)
    print("%f" % len(cumulative_df))
    cum_after = cumulative_df[columns].iloc[after_peak_day]
    cum_rate_after = cum_after/new_sizes
    active_rate_after = I[columns].iloc[after_peak_day]/new_sizes
    
    df  = pd.DataFrame({'model_tag': base_dir,
                        'cumulative_infected_40_days': np.round(cum_total_40_days), 
                        'cumulative_rate_40_days': cum_rate_40_days,
                        'active_rate_before_peak': active_rate_before,
                        'cumulative_before_peak': np.round(cum_before),
                        'cumulative_rate_before_peak': cum_rate_before,
                        'days_peak': days_peak,
                        'peak_active_infected_rate': peak_infected.values/new_sizes,                       
                        'cumulative_peak': np.round(cum_at_peak),
                        'cumulative_rate_peak': cum_at_peak_rate,
                        'active_rate_after_peak': active_rate_after,
                        'cumulative_after_peak': cum_after,
                        'cumulative_rate_after_peak': cum_rate_after,
                        'cumulative_infected_120_days': np.round(cum_total_100_days),
                        'cumulative_rate_120_days': np.round(cum_total_100_days/new_sizes)
                     })
    df['model_name'] = name
    active_df = I[columns].add_prefix('active_')
    cumulative_df = cumulative_df.add_prefix('cumulative_')
    return df
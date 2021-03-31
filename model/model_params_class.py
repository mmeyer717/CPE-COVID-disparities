"""Define class ModelParams which has relevant attributes

Params
------
transmission_rate: float
sip_start_date: int
initial_infection_multiplier: float
prison_infection_rate = None: float
police_contact_rate = None: float
police_group_size = int

Returns
-------
an object of the class ModelParams that has attributes for each of the above
parameters
"""
class Policy2Params:
    def __init__(self,
                 prison_sip_i_white, prison_sip_i_black,
                 jail_sip_i_white, jail_sip_i_black):
#         self.jail_release_shrink = jail_release_shrink
#         self.jail_release_date = jail_release_date,
        self.prison_sip_i_white = prison_sip_i_white
        self.prison_sip_i_black = prison_sip_i_black
        self.jail_sip_i_black = jail_sip_i_black
        self.jail_sip_i_white = jail_sip_i_white
        
class ModelParams:
    def __init__(self,transmission_rate, sip_start_date,initial_infection_multiplier,
                 prison_peak_date, prison_infection_rate,
                police_contact_rate = None,
                 police_group_size = None,
                 post_sip_t = None):
        self.transmission_rate = transmission_rate
        if post_sip_t == None:
            self.post_sip_transmission_rate = transmission_rate
        else:
            self.post_sip_transmission_rate = post_sip_t
            
        self.sip_start_date = sip_start_date
        self.prison_peak_date = prison_peak_date
        self.initial_infection_multiplier = initial_infection_multiplier
        self.prison_infection_rate = prison_infection_rate
        self.police_contact_rate = police_contact_rate
        self.police_group_size = police_group_size
        
    """Add attributes based on parameters we modified for uncertainty calculations
    """    
    def add_uncertainty_params(self, infection, pc, pgrp):
        self.prison_infection_rate = infection
        self.police_contact_rate = pc
        self.police_group_size = pgrp
        return self
        
    """Format name that is suitable for these model params
    """
    def get_name(self):
        suffix = '_pc%s_pgrp%s'%(self.police_contact_rate, self.police_group_size)
        name = 'I%s_PI%s_L%s'%(self.initial_infection_multiplier,
                               self.prison_infection_rate, self.sip_start_date) + \
                                   "__" + suffix
        return name
        
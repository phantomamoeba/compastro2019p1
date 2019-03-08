"""
Photon class ... handles propogation, pathing, absorption, escape, etc of a photon
"""

import utilities

class Photon:
    def __init__(self,wavelength,seed=None):
        self.w = wavelength
        self.seed = seed
        self.status = 0 # -1 = absorbed, 0 = propogating, 1 = free (escaped)




    def is_abosorbed(self):
        """
        Perform a random draw to see if this photon is aborbed at this event

        :return:
        """
        p = utilities.get_uniform_prob(seed=self.seed)

        if p < utilities.prob_absorb(self.w):
            self.status = -1
            return True
        else:
            return False


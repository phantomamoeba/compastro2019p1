"""
The simulation
"""

import photon
import utilities
import numpy as np
import matplotlib.pyplot as plt
import bisect

h = 6.626e-27
k = 1.381e-16
T = 50000
c = 2.998e10

class Blackbody_CDF:
    def __init__(self):
        # Read blackbody cumulative dist function (in terms of unitless energy term hv/kT)
        cdf, x_vals = [], []
        with open("cdf.txt") as f:
            x = 0.
            increment = 0.01
            for line in f:
                x += increment
                x = round(x,2)
                cdf.append(float(line.strip()))
                x_vals.append(x)
                if x==10.:
                    increment = 0.1
        self.x_vals = np.array(x_vals)
        self.cdf = np.array(cdf) / cdf[-1]
        np.random.seed(32)

    def gen_lambda(self):
        'Generate wavelengths in microns'
        # rnd drawn from uniform distribution
        rnd = np.random.random()
        cdf = self.cdf
        x_vals = self.x_vals
        # find location of rnd within the blackbody cdf
        i = bisect.bisect_right(cdf, rnd)
        # interpolate corresponding wavelength
        x = x_vals[-1]
        if i:
            # t is the interpolation amount
            t = (rnd - cdf[i-1]) / (cdf[i] - cdf[i-1])
            # interpolate linearly
            x = x_vals[i-1] + t * (x_vals[i]-x_vals[i-1])
        # convert unitless energy to wavelength in microns
        return h*c / (x*k*T) * 1e4

    def gen_spectrum(self, num_photons, cutoff_microns):
        'Return spectrum in microns'
        lambdas = []
        for _ in range(num_photons):
            l = self.gen_lambda()
            if l < cutoff_microns:
                lambdas.append(l)
        return lambdas

def main():
    #uncomment to set seed and make repeateable
    #utilities.prand_seed(137)

    blackbody = Blackbody_CDF()
    wavelengths = blackbody.gen_spectrum(num_photons=1000, cutoff_microns=0.5)

    print("Creating ....")
    photons = []
    for w in wavelengths:
        photons.append(photon.Photon(w))
    photons = np.array(photons)

    print("Starting propagation ...")
    escaped_wavelengths = [] 
    for i in range(len(wavelengths)):
        #todo: optionally re-seed
        utilities.prand_seed(None) #to use time

        while photons[i].status == 0:
            photons[i].propagate()
        if photons[i].status == 1:
            escaped_wavelengths.append(photons[i].w)

    plt.close('all')
    bins = np.linspace(0, 0.5, 100)
    plt.hist(wavelengths, bins, alpha=0.5)
    plt.hist(escaped_wavelengths, bins, alpha=0.5)
    plt.show()


if __name__ == '__main__':
    main()


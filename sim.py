"""
The simulation
"""

import photon
import utilities
import numpy as np
import matplotlib.pyplot as plt

def main():

    #uncomment to set seed and make repeateable
    #utilities.prand_seed(137)

    #wavelengths = np.linspace(0.01,3.0,1000) #10AA-30,000AA
    wavelengths = np.logspace(-2, 0.5, 100)  # 10AA ~ 30,000AA

    #photons = [photon.Photon(w) for w in wavelengths]*100

    #5000 wavelengths, 100 photons each
    #photons = np.array([[]*100 for w in wavelengths])

    print("Creating ....")
    photons = []
    for w in wavelengths:
        same_w = []
        for i in range(100):
            same_w.append(photon.Photon(w))
        photons.append(same_w)

    photons = np.array(photons)

    print("Starting propogation ...")

    #build up f_esc as we go
    f_esc = np.zeros(len(wavelengths))
    for i in range(len(wavelengths)):
        ct = 0
        for j in range(len(photons[i])):
            #todo: optionally re-seed
            utilities.prand_seed(None) #to use time

            while photons[i][j].status == 0:
                photons[i][j].propagate()
            if photons[i][j].status == 1:
                ct += 1

        f_esc[i] = float(ct)/len(photons[i])

        print("wavelength: %f [microns] f_esc = %0.1f" % (wavelengths[i],f_esc[i]*100.))



    #simple sanity check on total escape
    # result = [p.status for p in photons.flatten()]
    # n,b,_ = plt.hist(result,bins=2)
    # print(n,b)
    # plt.show()

    #results by wavelength
    # f_esc = np.zeros(len(wavelengths))
    # for i in range(len(f_esc)):
    #     f_esc[i] = [p.status for p in photons[i]].count(1)/len(photons[i])

    plt.close('all')
    plt.plot(wavelengths,f_esc)
    plt.xscale('log')
    plt.show()


if __name__ == '__main__':
    main()


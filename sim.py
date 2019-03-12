"""
The simulation
"""

import photon
import utilities
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

def plot(fn=None):
    #just plot
    wavelengths = np.load("wavelengths.npy")

    simruns = glob.glob("simrun_*.npy")
    sims = []
    for s in simruns:
        try:
            f_esc = np.load(s)
            sims.append(f_esc)

            run_id = int(s.split("_")[1].split('.'))
            if run_id > next_run_id:
                next_run_id = run_id
        except:
            pass

    f_esc = np.zeros(len(wavelengths))
    f_esc_err = np.zeros(len(wavelengths))

    sims = np.array(sims)

    for i in range(len(f_esc)):
        bin = sims[:, i]
        f_esc[i] = np.mean(bin)
        f_esc_err[i] = np.std(bin)

    plt.close('all')
    plt.plot(wavelengths,f_esc,label="mean f_esc")
    plt.fill_between(wavelengths,f_esc-f_esc_err,f_esc+f_esc_err,color='k',alpha=0.3,label=r"1-$\sigma$")
    plt.xscale('log')
    plt.legend()
    plt.title("f_esc by wavelength (%d simulations)" %(len(simruns)))
    plt.xlabel("wavelength bin [microns]")
    plt.ylabel("fraction of escaped photons")

    if fn is None:
        plt.show()
    else:
        plt.savefig(fn)

def progress_bar(percent, barLen = 20):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.2f%%" % (progress, percent * 100))
    sys.stdout.flush()

def main():

    #uncomment to set seed and make repeateable
    #utilities.prand_seed(137)

    #wavelengths = np.linspace(0.01,3.0,1000) #10AA-30,000AA

    #load the wavelenghts if they exist, otherwise create and save them (for future runs)
    try:
        wavelengths = np.load("wavelengths.npy")
    except:
        wavelengths = np.logspace(-2, 0.5, 1000)  # 10AA ~ 30,000AA
        np.save("wavelengths",wavelengths)

    len_waves = len(wavelengths)

    #photons = [photon.Photon(w) for w in wavelengths]*100

    #5000 wavelengths, 100 photons each
    #photons = np.array([[]*100 for w in wavelengths])


    #load all previous runs
    next_run_id = 0
    simruns = glob.glob("simrun_*.npy")
    sims = []
    for s in simruns:
        try:
            f_esc = np.load(s)
            sims.append(f_esc)

            run_id = int(s.split("_")[1].split('.')[0])
            if run_id > next_run_id:
                next_run_id = run_id
        except:
            pass

    next_run_id += 1

    print("Creating 100x photon packets [%f, %f microns] in %d bins ..." %(wavelengths[0], wavelengths[-1],len_waves))
    photons = []
    for w in wavelengths:
        same_w = []
        for i in range(100):
            same_w.append(photon.Photon(w))
        photons.append(same_w)

    photons = np.array(photons)

    print("Propagating simulation #%d ..." %(next_run_id) )

    #build up f_esc as we go
    f_esc = np.zeros(len_waves)
    for i in range(len_waves):
        ct = 0
        for j in range(len(photons[i])):
            #todo: optionally re-seed
            utilities.prand_seed(None) #to use time

            while photons[i][j].status == 0:
                photons[i][j].propagate()
            if photons[i][j].status == 1:
                ct += 1

        f_esc[i] = float(ct)/len(photons[i])

        #print("wavelength: %f [microns] f_esc = %0.1f" % (wavelengths[i],f_esc[i]*100.))
        progress_bar(float(i)/len_waves)

    progress_bar(1.0)
    print("\n")
    sims.append(f_esc)
    np.save("simrun_%d" %next_run_id, f_esc)

    #rebuild f_esc from all
    f_esc = np.zeros(len_waves)
    f_esc_err = np.zeros(len_waves)

    sims = np.array(sims)
    for i in range(len(f_esc)):
        bin = sims[:,i]
        f_esc[i] = np.mean(bin)
        f_esc_err[i] = np.std(bin)




    #simple sanity check on total escape
    # result = [p.status for p in photons.flatten()]
    # n,b,_ = plt.hist(result,bins=2)
    # print(n,b)
    # plt.show()

    #results by wavelength
    # f_esc = np.zeros(len(wavelengths))
    # for i in range(len(f_esc)):
    #     f_esc[i] = [p.status for p in photons[i]].count(1)/len(photons[i])

    # plt.close('all')
    # plt.plot(wavelengths,f_esc)
    # plt.fill_between(wavelengths,f_esc-f_esc_err,f_esc+f_esc_err,color='k',alpha=0.3)
    # plt.xscale('log')
    # plt.show()


if __name__ == '__main__':
    main()


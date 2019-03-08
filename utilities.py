"""
Collection of utilities such as random selection of direction, conversion between coord systems, validations plotting
"""
import numpy as np
import matplotlib.pyplot as plt

def get_direction(seed=None):
    """
    return a single coord tuple (theta, phi) where theta 0-pi and phi 0-2pi
    :param seed: if provided, will (re) seed the pRNG
    :return:
    """

    if seed is not None:
        np.random.seed(seed)

    #generate random uniform -1 to 1 as mu = cos(theta) ... so theta == arccos(mu)
    theta = np.arccos(2. * np.random.random() -1.)
    #phi runs 0 - 2pi
    phi = 2. * np.pi * np.random.random()

    return theta,phi



def test_direction(count=10000,bins=50,fn=None):
    """
    Plot histogram of directions to show uniform distribution over theta (as cos(mu)) and phi
    :param count: number of directions to get
    :param bins: number of bins
    :param fn: savefig to file, if provided
    :return:
    """
    theta = []
    phi = []

    for i in range(count):
        t,p = get_direction()
        theta.append(t)
        phi.append(p)

    theta = np.array(theta)
    phi = np.array(phi)

    plt.figure(figsize=(10, 3))
    plt.subplots_adjust(wspace=0.4, hspace=0.4)


    plt.subplot(131)
    plt.suptitle("")
    plt.title(r"$\theta$")
    plt.ylabel("count")
    plt.xlabel("radians")
    plt.hist(theta,bins=bins)


    plt.subplot(132)
    plt.suptitle("")
    plt.title(r"$\mu$ = cos($\theta$)")
    plt.ylabel("count")
    plt.xlabel("radians")
    plt.hist(np.cos(theta),bins=bins)

    plt.subplot(133)
    plt.suptitle("")
    plt.title(r"$\phi$")
    plt.ylabel("count")
    plt.xlabel("radians")
    plt.hist(phi,bins=bins)

    if fn is None:
        plt.show()
    else:
        plt.savefig(fn)


def sphere2cart(r,t,p):
    """
    Convert spherical polar coords to Cartessian
    :param r: radius
    :param t: polar angle (0-pi)
    :param p: azimuthal angle (0-2pi)
    :return: x,y,z
    """

    x = r * np.sin(t) * np.cos(p)
    y = r * np.sin(t) * np.sin(p)
    z = r * np.cos(t)

    return x,y,z



def prob_absorb(w):
    """

    Goes to constant ~ 20% for wavelengths < 0.05micron.
    Approaches 0% on longer wavelegnths (0.06% by 1 micron).

    :param w: wavelength in microns
    :param density:
    :return:
    """

    density_H = 1e-19 #g cm^-3
    density_dust = 1e-3 * density_H
    k_Th = (0.67e-24)/(1.67e-24) #cm^2/g Thomson opacity
    k_0 = 100. #cm^2 g^-1
    w_0 = 0.05 #microns
    k_d = k_0 * (w/w_0)**(-2.) if (w > w_0) else k_0

    return density_dust*k_d / (density_dust * k_d + density_H * k_Th)


def test_prob_absorb(fn=None):
    """
    Show the prob of absorption curve.

    :param fn: savefig to file, if provided
    :return:
    """

    p = []
    w = np.logspace(-3, 2, 10000)
    for x in w:
        p.append(prob_absorb(x))

    plt.title("Probability (0-1) of Absorption by Wavelength")
    plt.plot(w,p)
    plt.xscale("log")
    plt.ylabel("Prob of Absorption")
    plt.xlabel(r"$\lambda$ [microns")

    if fn is None:
        plt.show()
    else:
        plt.savefig(fn)
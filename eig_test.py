import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt

def sklearn_kde(data, x, bandwidth):
    kernel = KernelDensity(bandwidth=bandwidth)
    kernel.fit(data[:, np.newaxis])
    pdf = np.exp(kernel.score_samples(x[:,np.newaxis]))
    return pdf
    
    

dat= np.loadtxt("bcl_test/eigenval.txt")
print dat.shape

ev =  dat[:,1]
ev_num = dat[:,0]



sub = ev > 0
ev_num = ev_num[sub]
ev     = ev[sub]
ev = np.sqrt(ev) * 5.308 # convert from 1/ps to cm-1

sub = ev < 5000
ev_num = ev_num[sub]
ev     = ev[sub]

ev_search_val_cm = 1634
ev_found_ind = np.abs(ev - ev_search_val_cm).argmin()
print "Found {}cm-1 ~ {} at index {}".format(
        ev[ev_found_ind], ev_search_val_cm, ev_num[ev_found_ind])

w = np.linspace(0, 4000, 1000)
DoS = sklearn_kde(ev, w, 10.)
plt.plot(w, DoS)
    

#hist, bins = np.histogram(ev, bins=200, density=True)
#plt.plot((bins[1:] + bins[:-1])/2., hist)
#plt.title("Density of states")
#plt.xlabel("hbar*w, cm-1")
#plt.ylabel("density")
plt.show()

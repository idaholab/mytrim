import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("vacancy_count_comparison.dat", skiprows = 2, delimiter = ",")


gs = gridspec.GridSpec(2, 4)
gs.update(wspace=0.5)
ax1 = plt.subplot(gs[0, :2], )
ax2 = plt.subplot(gs[0, 2:])
ax3 = plt.subplot(gs[1, 1:3])

ax1.set_title("Si on C")
ax1.loglog(data[:,0], data[:,1], "b-", label="mytrim, mKP")
ax1.loglog(data[:,0], data[:,2], "g-",label="mytrim, detailed")
ax1.loglog(data[:,0], data[:,3], "r-", label="trim, mKP")
ax1.loglog(data[:,0], data[:,4], "m-",label="trim, detailed")
ax1.set_xlabel("Energy (keV)")
ax1.set_ylabel("Vacancies per ion")
ax1.legend(loc=0)

ax2.set_title("Xe on U")
ax2.loglog(data[:,0], data[:,5], "b-", label="mytrim, mKP")
ax2.loglog(data[:,0], data[:,6], "g-",label="mytrim, detailed")
ax2.loglog(data[:,0], data[:,7], "r-", label="trim, mKP")
ax2.loglog(data[:,0], data[:,8], "m-",label="trim, detailed")
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel("Vacancies per ion")

ax3.set_title("Cu on Cu")
ax3.loglog(data[:,0], data[:,9], "b-", label="mytrim, mKP")
ax3.loglog(data[:,0], data[:,10], "g-",label="mytrim, detailed")
ax3.loglog(data[:,0], data[:,11], "r-", label="trim, mKP")
ax3.loglog(data[:,0], data[:,12], "m-",label="trim, detailed")
ax3.set_xlabel("Energy (keV)")
ax3.set_ylabel("Vacancies per ion")

plt.show()

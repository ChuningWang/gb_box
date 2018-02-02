import netCDF4 as nc
import matplotlib.pyplot as plt

fh = nc.Dataset('./data/inverse.nc', 'r')
mt = fh.variables['t'][:]
s = fh.variables['s_box'][:]
f = fh.variables['f_box'][:]
sp = fh.variables['sp'][:]
Q02 = fh.variables['Q02'][:]
Q12 = fh.variables['Q12'][:]
Q23 = fh.variables['Q23'][:]
Q3p = fh.variables['Q3p'][:]
Wui = fh.variables['Wui'][:]
Wid = fh.variables['Wid'][:]

ax = plt.subplot(411)
plt.plot(mt, s[0, 0, :], '-', lw=2)
plt.plot(mt, s[1, 0, :], '-', lw=2)
plt.plot(mt, s[2, 0, :], '-', lw=2)
plt.plot([-10, -5], [0, 0], '--', lw=2)
plt.legend(('Upper', 'Intermediate', 'Deep', 'Pacific'), fontsize=6, loc=3)
plt.text(330, 28, 'BOX0')
plt.xlim(0, 366)
plt.ylim(27, 33)
plt.ylabel('Salinity [PSU]', fontsize=15)
ax.tick_params(axis='both', labelsize=10)

ax2 = ax.twinx()
ax2.plot(mt, f[0, :], 'grey')
plt.xlim(0, 366)
plt.ylabel(r'River Discharge [m$^3$s$^{-1}$]    ', fontsize=15)
ax2.spines['right'].set_color('grey')
ax2.yaxis.label.set_color('grey')
ax2.tick_params(axis='y', colors='grey')
ax2.tick_params(axis='both', labelsize=8)


ax = plt.subplot(412)
plt.plot(mt, s[0, 1, :], '-', lw=2)
plt.plot(mt, s[1, 1, :], '-', lw=2)
plt.plot(mt, s[2, 1, :], '-', lw=2)
plt.text(330, 28, 'BOX1')
plt.xlim(0, 366)
plt.ylim(27, 33)
ax.tick_params(axis='both', labelsize=8)

ax2 = ax.twinx()
ax2.plot(mt, f[1, :], 'grey')
plt.xlim(0, 366)
ax2.spines['right'].set_color('grey')
ax2.yaxis.label.set_color('grey')
ax2.tick_params(axis='y', colors='grey')
ax2.tick_params(axis='both', labelsize=8)


ax = plt.subplot(413)
plt.plot(mt, s[0, 2, :], '-', lw=2)
plt.plot(mt, s[1, 2, :], '-', lw=2)
plt.plot(mt, s[2, 2, :], '-', lw=2)
plt.text(330, 28, 'BOX2')
plt.xlim(0, 366)
plt.ylim(27, 33)
ax.tick_params(axis='both', labelsize=8)

ax2 = ax.twinx()
ax2.plot(mt, f[2, :], 'grey')
plt.xlim(0, 366)
ax2.spines['right'].set_color('grey')
ax2.yaxis.label.set_color('grey')
ax2.tick_params(axis='y', colors='grey')
ax2.tick_params(axis='both', labelsize=8)


ax = plt.subplot(414)
plt.plot(mt, s[0, 3, :], '-', lw=2)
plt.plot(mt, s[1, 3, :], '-', lw=2)
plt.plot(mt, s[2, 3, :], '-', lw=2)
plt.plot(mt, sp, '--', lw=2)
plt.text(330, 28, 'BOX3')
plt.xlim(0, 366)
plt.ylim(27, 33)
plt.xlabel('Yearday')
ax.tick_params(axis='both', labelsize=8)

ax2 = ax.twinx()
ax2.plot(mt, f[3, :], 'grey')
plt.xlim(0, 366)
ax2.spines['right'].set_color('grey')
ax2.yaxis.label.set_color('grey')
ax2.tick_params(axis='y', colors='grey')
ax2.tick_params(axis='both', labelsize=8)


plt.savefig('s_box.eps', format='eps')
plt.close()

# --------------------------------------------------------------------------------------

# Smooth inverse model outputs
from scipy.signal import butter, filtfilt
fs = 1./(mt[1]-mt[0])  # Sampling frequency [day^-1]
cutoff = 1./30  # Cutoff frequency [day^-1]
b, a = butter(5, cutoff/(0.5*fs))

ds0 = filtfilt(b, a, s[0, 2, :]-s[0, 0, :])
ds1 = filtfilt(b, a, s[0, 2, :]-s[0, 1, :])
ds2 = filtfilt(b, a, s[0, 3, :]-s[0, 2, :])
ds3 = filtfilt(b, a, sp-s[0, 3, :])
dsi0 = filtfilt(b, a, s[1, 0, :]-s[0, 0, :])
dsi1 = filtfilt(b, a, s[1, 1, :]-s[0, 1, :])
dsi2 = filtfilt(b, a, s[1, 2, :]-s[0, 2, :])
dsi3 = filtfilt(b, a, s[1, 3, :]-s[0, 3, :])
dsd0 = filtfilt(b, a, s[2, 0, :]-s[1, 0, :])
dsd1 = filtfilt(b, a, s[2, 1, :]-s[1, 1, :])
dsd2 = filtfilt(b, a, s[2, 2, :]-s[1, 2, :])
dsd3 = filtfilt(b, a, s[2, 3, :]-s[1, 3, :])
# dsd = filtfilt(b, a, inv['sd']-inv['si0'])

# ---------------------------------------------------------------

fig, axarr = plt.subplots(4, sharex=True)

plt.sca(axarr[0])
plt.plot(mt, Q02, 'k', lw=2)
plt.ylabel(r'Q$_{02}$ [m$^3$s$^{-1}$]')

plt.sca(axarr[1])
plt.plot(mt, Q12, 'k', lw=2)
plt.ylabel(r'Q$_{12}$')

plt.sca(axarr[2])
plt.plot(mt, Q23, 'k', lw=2)
plt.ylabel(r'Q$_{23}$')

plt.sca(axarr[3])
plt.plot(mt, Q3p, 'k', lw=2)
plt.ylabel(r'Q$_{3p}$')
plt.xlabel('Yearday')
plt.xlim(0, 366)

pltds = 0
if pltds == 1:
    for i in range(4):
        ax2 = axarr[i].twinx()
        if i == 0:
            ax2.plot(mt, ds0, 'r')
            plt.ylabel(r'S$_{2u}$-S$_{0u}$ [PSU]')
        elif i == 1:
            ax2.plot(mt, ds1, 'r')
            plt.ylabel(r'S$_{2u}$-S$_{1u}$')
        elif i == 2:
            ax2.plot(mt, ds2, 'r')
            plt.ylabel(r'S$_{3u}$-S$_{2u}$')
        elif i == 3:
            ax2.plot(mt, ds3, 'r')
            plt.ylabel(r'S$_{p}$-S$_{3u}$')
        plt.xlim(0, 366)
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')

pltf = 1
if pltf == 1:
    for i in range(4):
        ax2 = axarr[i].twinx()
        ax2.plot(mt, f[i, :], 'r')
        plt.xlim(0, 366)
        if i == 0:
            plt.ylabel(r'F [m$^3$s$^{-1}$]')
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')




# plt.show(block=False)
plt.savefig('Q3_f.eps', format='eps')
plt.close()


# ---------------------------------------------------------------

fig, axarr = plt.subplots(4, sharex=True)

plt.sca(axarr[0])
plt.plot(mt, Wui[0], 'k', lw=2)
plt.ylabel(r'W$_{ui}$ [m$^3$s$^{-1}$]')

plt.sca(axarr[1])
plt.plot(mt, Wui[1], 'k', lw=2)

plt.sca(axarr[2])
plt.plot(mt, Wui[2], 'k', lw=2)

plt.sca(axarr[3])
plt.plot(mt, Wui[3], 'k', lw=2)
plt.xlabel('Yearday')
plt.xlim(0, 366)

pltds = 1
if pltds == 1:
    for i in range(4):
        ax2 = axarr[i].twinx()
        if i == 0:
            ax2.plot(mt, dsi0, 'r')
            plt.ylabel(r'S$_{0i}$-S$_{0u}$ [PSU]')
        elif i == 1:
            ax2.plot(mt, dsi1, 'r')
            plt.ylabel(r'S$_{1i}$-S$_{1u}$')
        elif i == 2:
            ax2.plot(mt, dsi2, 'r')
            plt.ylabel(r'S$_{2i}$-S$_{2u}$')
        elif i == 3:
            ax2.plot(mt, dsi3, 'r')
            plt.ylabel(r'S$_{3i}$-S$_{3u}$')
        plt.xlim(0, 366)
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')


plt.savefig('Wui_s.eps', format='eps')
plt.close()

# ---------------------------------------------------------------

fig, axarr = plt.subplots(4, sharex=True)

plt.sca(axarr[0])
plt.plot(mt, Wid[0], 'k', lw=2)
plt.ylabel(r'W$_{id}$ [m$^3$s$^{-1}$]')

plt.sca(axarr[1])
plt.plot(mt, Wid[1], 'k', lw=2)

plt.sca(axarr[2])
plt.plot(mt, Wid[2], 'k', lw=2)

plt.sca(axarr[3])
plt.plot(mt, Wid[3], 'k', lw=2)
plt.xlabel('Yearday')
plt.xlim(0, 366)

pltds = 1
if pltds == 1:
    for i in range(4):
        ax2 = axarr[i].twinx()
        if i == 0:
            ax2.plot(mt, dsd0, 'r')
            plt.ylabel(r'S$_{0i}$-S$_{0u}$ [PSU]')
        elif i == 1:
            ax2.plot(mt, dsd1, 'r')
            plt.ylabel(r'S$_{1i}$-S$_{1u}$')
        elif i == 2:
            ax2.plot(mt, dsd2, 'r')
            plt.ylabel(r'S$_{2i}$-S$_{2u}$')
        elif i == 3:
            ax2.plot(mt, dsd3, 'r')
            plt.ylabel(r'S$_{3i}$-S$_{3u}$')
        plt.xlim(0, 366)
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')


plt.savefig('Wid_s.eps', format='eps')
plt.close()

# ----------------------------------------------------------------------------------------------
import numpy as np

cal_inverse = 1
if cal_inverse == 1:
    from box_gb import box_inverse
    inv = box_inverse.box_inverse(cal_clim=1, svpth='./data/')

mt = inv['mt']
s = inv['s_box']
f = inv['f_box']
sp = inv['sp']
Q02 = inv['Q02']
Q12 = inv['Q12']
Q23 = inv['Q23']
Q3p = inv['Q3p']
Wui0 = inv['Wui0']
Wui3 = inv['Wui3']
Wid = inv['Wid']

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
dsd = filtfilt(b, a, inv['sd']-inv['si0'])

mt2 = 0.5*(mt[:-1]+mt[1:])

import matplotlib.pyplot as plt

ax = plt.subplot(411)
plt.plot(mt, Q02, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 0, :]), 'r')
ax2.plot(mt, ds0, 'r')
plt.xlim(0, 366)

ax = plt.subplot(412)
plt.plot(mt, Q12, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 1, :]), 'r')
ax2.plot(mt, ds1, 'r')
plt.xlim(0, 366)

ax = plt.subplot(413)
plt.plot(mt, Q23, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 3, :]-s[0, 2, :]), 'r')
ax2.plot(mt, ds2, 'r')
plt.xlim(0, 366)

ax = plt.subplot(414)
plt.plot(mt, Q3p, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(sp[:]-s[0, 3, :]), 'r')
ax2.plot(mt, ds3, 'r')
plt.xlim(0, 366)


# plt.show(block=False)
plt.savefig('Q3.eps', format='eps')
plt.close()

# -----------------------------------------------------------------------------------------

dds0 = np.diff(ds0)/np.diff(mt)
dds1 = np.diff(ds1)/np.diff(mt)
dds2 = np.diff(ds2)/np.diff(mt)
dds3 = np.diff(ds3)/np.diff(mt)


ax = plt.subplot(411)
plt.plot(mt, Q02, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 0, :]), 'r')
ax2.plot(mt[:-1], dds0, 'r')
plt.xlim(0, 366)

ax = plt.subplot(412)
plt.plot(mt, Q12, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 1, :]), 'r')
ax2.plot(mt[:-1], dds1, 'r')
plt.xlim(0, 366)

ax = plt.subplot(413)
plt.plot(mt, Q23, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 3, :]-s[0, 2, :]), 'r')
ax2.plot(mt[:-1], dds2, 'r')
plt.xlim(0, 366)

ax = plt.subplot(414)
plt.plot(mt, Q3p, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(sp[:]-s[0, 3, :]), 'r')
ax2.plot(mt[:-1], dds3, 'r')
plt.xlim(0, 366)


# plt.show(block=False)
plt.savefig('Q4.eps', format='eps')
plt.close()

# -----------------------------------------------------------------------------------------

ax = plt.subplot(411)
plt.plot(mt, Q02, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 0, :]), 'r')
ax2.plot(mt, f[0, :], 'r')
plt.xlim(0, 366)

ax = plt.subplot(412)
plt.plot(mt, Q12, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 1, :]), 'r')
ax2.plot(mt, f[1, :], 'r')
plt.xlim(0, 366)

ax = plt.subplot(413)
plt.plot(mt, Q23, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(s[0, 3, :]-s[0, 2, :]), 'r')
ax2.plot(mt, f[2, :], 'r')
plt.xlim(0, 366)

ax = plt.subplot(414)
plt.plot(mt, Q3p, '--', lw=2)

ax2 = ax.twinx()
# ax2.plot(mt2, np.diff(sp[:]-s[0, 3, :]), 'r')
ax2.plot(mt, f[3, :], 'r')
plt.xlim(0, 366)


# plt.show(block=False)
plt.savefig('Q5.eps', format='eps')
plt.close()

# -----------------------------------------------------------------------------------------




plt.figure()
ax = plt.subplot(311)
plt.plot(mt, Wui0)

ax2 = ax.twinx()
ax2.plot(mt, dsi0, 'r')
ax2.plot(mt, dsi1, 'k')
ax2.plot(mt, dsi2, 'g')

ax = plt.subplot(312)
plt.plot(mt, Wui3)

ax2 = ax.twinx()
ax2.plot(mt, dsi3, 'r')

ax = plt.subplot(313)
plt.plot(mt, Wid)

ax2 = ax.twinx()
ax2.plot(mt, dsd, 'r')

plt.savefig('W1.eps', format='eps')
plt.close()

ax = plt.subplot(411)
plt.plot(mt, s[0, 0, :], '-', lw=2)
plt.plot(mt, s[1, 0, :], '-', lw=2)
plt.plot(mt, s[2, 0, :], '-', lw=2)
plt.xlim(0, 366)
plt.ylim(28, 33)
plt.legend('Upper', 'Intermediate', 'Deep')

ax2 = ax.twinx()
ax2.plot(mt, f[0, :], 'grey')
plt.xlim(0, 366)


plt.subplot(412)
plt.plot(mt, s[0, 1, :], '-', lw=2)
plt.plot(mt, s[1, 1, :], '-', lw=2)
plt.plot(mt, s[2, 1, :], '-', lw=2)
plt.xlim(0, 366)
plt.ylim(28, 33)

plt.subplot(413)
plt.plot(mt, s[0, 2, :], '-', lw=2)
plt.plot(mt, s[1, 2, :], '-', lw=2)
plt.plot(mt, s[2, 2, :], '-', lw=2)
plt.xlim(0, 366)
plt.ylim(28, 33)

plt.subplot(414)
plt.plot(mt, s[0, 3, :], '-', lw=2)
plt.plot(mt, s[1, 3, :], '-', lw=2)
plt.plot(mt, s[2, 3, :], '-', lw=2)
plt.plot(mt, sp, '--', lw=2)
plt.xlim(0, 366)
plt.ylim(28, 33)
plt.xlabel('Yearday')
plt.ylabel('Salinity [PSU]')

plt.savefig('s_box.eps', format='eps')
plt.close()


# import netCDF4 as nc
# f = nc.Dataset('./data/inverse.nc', 'r')
# 
# mt = f.variables['t'][:]
# Q02 = f.variables['Q02'][:]
# Q12 = f.variables['Q12'][:]
# Q23 = f.variables['Q23'][:]
# Q3p = f.variables['Q3p'][:]
# Wui = f.variables['Wui'][:]
# Wid = f.variables['Wid'][:]
# s = f.variables['s_box'][:]
# F = f.variables['f_box'][:]
# sp = f.variables['sp'][:]
# 
# f.close()
# 
# # Smooth inverse model outputs
# from scipy.signal import butter, filtfilt
# fs = 1./(mt[1]-mt[0])  # Sampling frequency [day^-1]
# cutoff = 1./30  # Cutoff frequency [day^-1]
# b, a = butter(5, cutoff/(0.5*fs))
# # Q02_f = filtfilt(b, a, Q02)
# # Q12_f = filtfilt(b, a, Q12)
# # Q23_f = filtfilt(b, a, Q23)
# # Q3p_f = filtfilt(b, a, Q3p)
# # Wui_f = filtfilt(b, a, Wui, axis=1)
# # Wid_f = filtfilt(b, a, Wid, axis=1)
# Q02_f = Q02
# Q12_f = Q12
# Q23_f = Q23
# Q3p_f = Q3p
# Wui_f = Wui
# Wid_f = Wid
# 
# ds0 = filtfilt(b, a, s[0, 2, :]-s[0, 0, :])
# ds1 = filtfilt(b, a, s[0, 2, :]-s[0, 1, :])
# ds2 = filtfilt(b, a, s[0, 3, :]-s[0, 2, :])
# ds3 = filtfilt(b, a, sp-s[0, 3, :])
# 
# mt2 = 0.5*(mt[:-1]+mt[1:])
# Q02_f2 = 0.5*(Q02_f[:-1]+Q02_f[1:])
# Q12_f2 = 0.5*(Q12_f[:-1]+Q12_f[1:])
# Q23_f2 = 0.5*(Q23_f[:-1]+Q23_f[1:])
# Q3p_f2 = 0.5*(Q3p_f[:-1]+Q3p_f[1:])
# dds0 = np.diff(ds0)/np.diff(mt)
# dds1 = np.diff(ds1)/np.diff(mt)
# dds2 = np.diff(ds2)/np.diff(mt)
# dds3 = np.diff(ds3)/np.diff(mt)

# ----------------------------------------------------------------------------------------------
# Make plots
pltfig = -1
if pltfig == 1:
    import matplotlib.pyplot as plt

    plt.subplot(411)
    plt.plot(mt, s[0, 0, :], '-', lw=2)
    plt.plot(mt, s[1, 0, :], '-', lw=2)
    plt.plot(mt, s[2, 0, :], '-', lw=2)
    plt.xlim(0, 366)
    plt.ylim(26, 32)

    plt.subplot(412)
    plt.plot(mt, s[0, 1, :], '-', lw=2)
    plt.plot(mt, s[1, 1, :], '-', lw=2)
    plt.plot(mt, s[2, 1, :], '-', lw=2)
    plt.xlim(0, 366)
    plt.ylim(26, 32)

    plt.subplot(413)
    plt.plot(mt, s[0, 2, :], '-', lw=2)
    plt.plot(mt, s[1, 2, :], '-', lw=2)
    plt.plot(mt, s[2, 2, :], '-', lw=2)
    plt.xlim(0, 366)
    plt.ylim(26, 32)

    plt.subplot(414)
    plt.plot(mt, s[0, 3, :], '-', lw=2)
    plt.plot(mt, s[1, 3, :], '-', lw=2)
    plt.plot(mt, s[2, 3, :], '-', lw=2)
    plt.plot(mt, sp, '--', lw=2)
    plt.xlim(0, 366)
    plt.ylim(30, 33)

    plt.savefig('s_box.eps', format='eps')
    plt.close()

    # ----------------------------------------------------------------------------------

    ax = plt.subplot(411)
    plt.plot(mt, Q02_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    # ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 0, :]), 'r')
    # ax2.plot(mt2, dds0, 'r')
    ax2.plot(mt, ds0, 'r')
    plt.xlim(0, 366)

    ax = plt.subplot(412)
    plt.plot(mt, Q12_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    # ax2.plot(mt2, np.diff(s[0, 2, :]-s[0, 1, :]), 'r')
    # ax2.plot(mt2, dds1, 'r')
    ax2.plot(mt, ds1, 'r')
    plt.xlim(0, 366)

    ax = plt.subplot(413)
    plt.plot(mt, Q23_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    # ax2.plot(mt2, np.diff(s[0, 3, :]-s[0, 2, :]), 'r')
    # ax2.plot(mt2, dds2, 'r')
    ax2.plot(mt, ds2, 'r')
    plt.xlim(0, 366)

    ax = plt.subplot(414)
    plt.plot(mt, Q3p_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    # ax2.plot(mt2, np.diff(sp[:]-s[0, 3, :]), 'r')
    # ax2.plot(mt2, dds3, 'r')
    ax2.plot(mt, ds3, 'r')
    plt.xlim(0, 366)


    # plt.show(block=False)
    plt.savefig('Q.eps', format='eps')
    plt.close()

    # ----------------------------------------------------------------------------------

    ax = plt.subplot(411)
    plt.plot(mt, Q02_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    ax2.plot(mt, F[0, :], 'r')
    plt.xlim(0, 366)

    ax = plt.subplot(412)
    plt.plot(mt, Q12_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    ax2.plot(mt, F[1, :], 'r')
    plt.xlim(0, 366)

    ax = plt.subplot(413)
    plt.plot(mt, Q23_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    ax2.plot(mt, F[2, :], 'r')
    plt.xlim(0, 366)

    ax = plt.subplot(414)
    plt.plot(mt, Q3p_f, '--', lw=2)
    # plt.ylim(0, 80000)

    ax2 = ax.twinx()
    ax2.plot(mt, F[3, :], 'r')
    plt.xlim(0, 366)


    # plt.show(block=False)
    plt.savefig('Q2.eps', format='eps')
    plt.close()

    # ----------------------------------------------------------------------------------

    plt.subplot(411)
    plt.plot(mt, Wui_f[0, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    plt.subplot(412)
    plt.plot(mt, Wui_f[1, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    plt.subplot(413)
    plt.plot(mt, Wui_f[2, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    plt.subplot(414)
    plt.plot(mt, Wui_f[3, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    # plt.show(block=False)
    plt.savefig('omega_ui.eps', format='eps')
    plt.close()

    plt.subplot(411)
    plt.plot(mt, Wid_f[0, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    plt.subplot(412)
    plt.plot(mt, Wid_f[1, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    plt.subplot(413)
    plt.plot(mt, Wid_f[2, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    plt.subplot(414)
    plt.plot(mt, Wid_f[3, :], '--', lw=2)
    plt.xlim(0, 366)
    # plt.ylim(0, 80000)

    # plt.show(block=False)
    plt.savefig('omega_id.eps', format='eps')
    plt.close()


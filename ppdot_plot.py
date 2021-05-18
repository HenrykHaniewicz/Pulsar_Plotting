#!/usr/local/bin/python3
# Plots a P-Pdot diagram (as PRESTO plotter no longer works for me!)

from builtins import zip
import numpy as np
import routines as ro
import pulsars
import matplotlib.pyplot as plt

# Use color?
usecolor = True

# Find a list of the "good" pulsars:  those not in GCs and with a measured pdot
# Also identify which pulsars are "special"
numgood = 0
numGC = 0
numpd0 = 0
ps = []
pds = []
rrats = []
radios = []
nonradios = []
magnetars = []
hepsrs = []
snrs = []
binaries = []
for psr in pulsars.psrs:
    # Ignore pulsars without measured Pdot
    if psr.pd==0.0:
        numpd0 += 1
        continue
    # Ignore globular cluster pulsars
    elif (psr.assoc is not None and 'GC' in psr.assoc):
        numGC += 1
        continue
    else:
        ps.append(psr.p)
        pds.append(psr.pd)
        if psr.type is not None:
            if 'RRAT' in psr.type: rrats.append(numgood)
            if 'NRAD' in psr.type:
                nonradios.append(numgood)
            if ('AXP' in psr.type) and ('NRAD' in psr.type):
                pass
            elif 'AXP' in psr.type:
                magnetars.append(numgood)
            if 'HE' in psr.type: hepsrs.append(numgood)
        if numgood not in nonradios:
            radios.append(numgood)
        if psr.assoc is not None:
            if 'SNR' in psr.assoc: snrs.append(numgood)
        if psr.binary:
            binaries.append(numgood)
        numgood += 1
ps = np.asarray(ps)
pds = np.asarray(pds)
rrats = np.asarray(rrats)
radios = np.asarray(radios)
nonradios = np.asarray(nonradios)
magnetars = np.asarray(magnetars)
hepsrs = np.asarray(hepsrs)
snrs = np.asarray(snrs)
binaries = np.asarray(binaries)

print( f"Plotting {numgood} pulsars total:" )
print( f"  {len(radios)} radio, {len(nonradios)} non-radio" )
print( f"  RRATs: {len(rrats)}" )
print( f"  magnetars: {len(magnetars)}" )
print( f"  high-energy: {len(hepsrs)}" )
print( f"  in SNRs: {len(snrs)}" )
print( f"  in binaries: {len(binaries)}" )
print( f"Rejected {numpd0} for having no p-dot and {numGC} for being in a cluster" )

# Now set up the plot
plims = np.asarray([0.001, 20.0])
pdlims = np.asarray([1e-22, 1e-9])
dpdpd = (np.log10(plims[1]) - np.log10(plims[0])) / (np.log10(pdlims[1]) - np.log10(pdlims[0]))
grey = '0.8'
greytext = '0.3'
plt.figure(num=None, figsize=(9, 9), dpi=200)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(plims)
ax.set_ylim(pdlims)
# Make period labels be non-scientific notation
ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter("%g"))

# Plot magnetic field lines
Bs_to_plot = [8, 10, 12, 14]
for logB in Bs_to_plot:
    plt.plot(plims, ro.pdot_from_B(plims, 10.0**logB), '-', color=grey)
    if logB==14:
        y = 4e-10
        x = ro.pdot_from_B(y, 10.0**logB)
    elif logB==8:
        x = 0.05
        y = ro.pdot_from_B(x, 10.0**logB)
    else:
        x = 1.1 * plims[0]
        y = ro.pdot_from_B(x, 10.0**logB)
    plt.text(x, y, "$10^{%d}$ G"%logB, color=greytext,
        horizontalalignment='left', verticalalignment='baseline',
        rotation=np.degrees(np.arctan(-1.0 * dpdpd)))

# Plot Edot lines
Edots_to_plot = [31, 34, 37, 40]
for logEdot in Edots_to_plot[::-1]:
    plt.plot(plims, ro.pdot_from_edot(plims, 10.0**logEdot), '-', color=grey)
    if logEdot > 31:
        y = 5e-10
        x = 0.6 * (y * 4e45 * np.pi * np.pi / 10.0**logEdot)**(1.0/3.0)
    else:
        y = 1e-21
        x = 0.6 * (y * 4e45 * np.pi * np.pi / 10.0**logEdot)**(1.0/3.0)
    plt.text(x, y, "$10^{%d}$ erg/s"%logEdot, color=greytext,
        horizontalalignment='left', verticalalignment='baseline',
        rotation=np.degrees(np.arctan(3.0 * dpdpd)))

# Plot Age lines
Ages_to_plot = [3, 5, 7, 9, 11]
Ages_labels = ['1 Kyr', '100 Kyr', '10 Myr', '1 Gyr', '100 Gyr']
for logAge, label in zip(Ages_to_plot, Ages_labels):
    plt.plot(plims, ro.pdot_from_age(plims, 10.0**logAge), '-', color=grey)
    x = 1.1 * plims[0]
    plt.text(x, 1.1 * ro.pdot_from_age(x, 10.0**logAge), label, color=greytext,
        horizontalalignment='left', verticalalignment='bottom',
        rotation=np.degrees(np.arctan(1.0 * dpdpd)))

# Now plot the radio pulsars as black dots
plt.plot(ps[radios], pds[radios], '.', color='0.3', ms=3, label = "Radio pulsars" )

# Plot the HE and non-radio pulsars as triagles
# Assume that all non-radio pulsars are high-energy emitters
all_he = np.unique(np.concatenate((hepsrs, nonradios)))
#color = 'magenta' if usecolor else 'black'
#plt.plot(ps[all_he], pds[all_he], '^', ms=6, mew=1.1, mec=color, mfc='none', label="X-ray/$\gamma$-ray")

# Plot the binaries as circles
plt.plot(ps[binaries], pds[binaries], 'ko', ms=8, mfc='none', label = "Binaries")


# J185*+00
plt.plot( 0.022838764007, 1.2011e-19, 'm*', label = "PSR J1851+0010", ms = 20 ) # J1851
plt.plot( 0.033402772370, 7.403e-19, 'r*', label = "PSR J1853+0008", ms = 20 ) # J1853
# J1936s
plt.plot( 0.05834509339955, 9.36e-20, 'g*', label = "PSR J1936+1805", ms = 20 ) # J1936+18
plt.plot( 0.0315888880147147, 4.6e-20, 'b*', label = "PSR J1936+2142", ms = 20 ) # J1936+21

# Plot the SNRs as stars
#color = 'darkorange' if usecolor else 'black'
#mew = 1.0 if usecolor else 0.7
#plt.plot(ps[snrs], pds[snrs], '*', ms=14, mfc='none', mew=mew, mec=color, label="SNR Assoc")

# Plot the magnetars as filled triangles
color = 'cyan' if usecolor else 'black'
plt.plot(ps[magnetars], pds[magnetars], '.', color='0.3', ms=3)

# Plot the RRATs as x's
color = 'red' if usecolor else 'black'
plt.plot(ps[rrats], pds[rrats], '.', color='0.3', ms=3)

plt.xlabel("Spin period (s)")
plt.ylabel(r"Spin period derivative (ss$^{-1}$)")

ax.legend(loc='lower right', numpoints=1)

plt.savefig( f"chapter4_ppdot_color_{pulsars.version}.pdf" if usecolor else f"chapter4_ppdot_{pulsars.version}.pdf", format = 'pdf')
#plt.show()

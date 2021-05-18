#!/usr/local/bin/python3
# Plots a P-Porb diagram (as PRESTO plotter no longer works for me!)

from builtins import zip
import numpy as np
import routines as ro
import pulsars
import matplotlib.pyplot as plt

# Use color?
usecolor = True

# Find a list of the "good" pulsars:  those not in GCs and with a measured pb
# Also identify which pulsars are "special"
numgood = 0
numGC = 0
numpb0 = 0
ps = []
pbs = []
rrats = []
radios = []
nonradios = []
magnetars = []
hepsrs = []
snrs = []
binaries = []
for psr in pulsars.psrs:
    # Ignore pulsars without measured PB
    try:
        tst = (psr.pb == None)
    except AttributeError:
        numpb0 += 1
        continue

    if psr.pb is 0.0:
        numpb0 += 1
        continue
    # Ignore globular cluster pulsars
    elif (psr.assoc is not None and 'GC' in psr.assoc):
        numGC += 1
        continue
    else:
        ps.append(psr.p)
        pbs.append(psr.pb)
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
pbs = np.asarray(pbs)
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
print( f"Rejected {numpb0} for having no orbital period and {numGC} for being in a cluster" )

# Now set up the plot
plims = np.asarray([0.01, 0.1])
pblims = np.asarray([0.01, 100])
dpbpb = (np.log10(plims[1]) - np.log10(plims[0])) / (np.log10(pblims[1]) - np.log10(pblims[0]))
grey = '0.8'
greytext = '0.3'
plt.figure(num=None, figsize=(9, 9), dpi=200)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(pblims)
ax.set_ylim(plims)
# Make period labels be non-scientific notation
ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter("%g"))


# Now plot the radio pulsars as black dots
plt.plot(pbs[radios], ps[radios], '.', color='0.3', ms=3, label = "Radio PSRs" )

# Plot the HE and non-radio pulsars as triagles
# Assume that all non-radio pulsars are high-energy emitters
all_he = np.unique(np.concatenate((hepsrs, nonradios)))
#color = 'magenta' if usecolor else 'black'
#plt.plot(ps[all_he], pbs[all_he], '^', ms=6, mew=1.1, mec=color, mfc='none', label="X-ray/$\gamma$-ray")

# Plot the binaries as circles
#plt.plot(pbs[binaries], ps[binaries], 'ko', ms=8, mfc='none', label = "Binaries")


# J185*+00
plt.plot( 2.003462165, 0.022838764007, 'm*', label = "PSR J1851+0010", ms = 20 ) # J1851
plt.plot( 9.6129413, 0.033402772370, 'y*', label = "PSR J1853+0008", ms = 20 ) # J1853
# J1936s
#plt.plot( 0.058345093399, 0.05834509339, 'go', label = "PSR J1936+1805" ) # J1936+18
plt.plot( 0.757, 0.0315888880147147, 'b*', label = "PSR J1936+2142", ms = 20 ) # J1936+21


pers = np.array([ 45.8, 76.5, 22.7, 62.5, 40.9, 37.9, 95.1, 28.5, 21.5, 104.2, 41.0, 27.3, 59.0, 185.5, 17.0 ])
pers = pers/1000.0

orbpers = [ 4.072, 0.38, 0.102, 2.616, 8.634, 0.421, 13.638, 0.32, 0.183, 18.779, 1.176, 0.206, 0.323, 45.060, 0.078 ]

plt.plot( orbpers, pers, '*', color = 'goldenrod', label = "DNSs"  )

# Plot the SNRs as stars
#color = 'darkorange' if usecolor else 'black'
#mew = 1.0 if usecolor else 0.7
#plt.plot(ps[snrs], pbs[snrs], '*', ms=14, mfc='none', mew=mew, mec=color, label="SNR Assoc")

# Plot the magnetars as filled triangles
#color = 'cyan' if usecolor else 'black'
#plt.plot(pbs[magnetars], ps[magnetars], '.', color='0.3', ms=3)

# Plot the RRATs as x's
#color = 'red' if usecolor else 'black'
#plt.plot(pbs[rrats], ps[rrats], '.', color='0.3', ms=3)

plt.xlabel("Orbital period (days)")
plt.ylabel(r"Spin period (s)")

ax.legend(loc='upper left', numpoints=1)

plt.savefig( f"chapter4_pporb_color_{pulsars.version}.pdf" if usecolor else f"chapter4_pporb_{pulsars.version}.pdf", format = 'pdf')
#plt.show()

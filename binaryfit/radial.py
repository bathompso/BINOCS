import sys, binocs
import numpy as np
import matplotlib.pyplot as plt

mbins = [[-1, 99], [0.2, 0.8], [0.8, 1.1], [1.1, 1.6], [1.6, 2.5], [2.5, 99]]

# Read options
options = binocs.readopt((sys.argv)[1])

# Determine cluster designation
dirsplit = options['data'].split('/')
namesplit = dirsplit[len(dirsplit)-1].split('.')
clnum = (namesplit[0])[3:len(namesplit[0])]

# Find cluster in Dias catalog
diaspath = "/Users/bathompso/Dropbox/dias.txt"
df = open(diaspath, 'r')
dias = df.read().splitlines()
df.close()
clidx = [x for x in range(len(dias)) if (dias[x])[0:19].find("NGC") >= 0 and (dias[x])[0:19].find(clnum) >= 0]
print "\nFound Cluster in Dias catalog at line", clidx[0]

# Convert Dias entries to decimal Ra/Dec
clcoo = [float(x) for x in (dias[clidx[0]])[18:38].split()]
clra = 15.0 * (clcoo[0] + clcoo[1]/60.0 + clcoo[2]/3600.0)
cldec = clcoo[3] + clcoo[4]/60.0 + clcoo[5]/3600.0
print "Dias cluster coordinates: %9.5f  %9.5f" % (clra, cldec)

# Read in results file
df = open(options['data']+'--results.txt', 'r')
lines = df.read().splitlines()
df.close()
mpri, msec, rad = np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines))
for l in range(len(lines)):
	tmp = lines[l].split()
	mpri[l] = float(tmp[22])
	msec[l] = float(tmp[24])
	rad[l] = np.sqrt( ((float(tmp[1]) - clra) * np.cos(cldec * np.pi/180.0))**2.0 + (float(tmp[2]) - cldec)**2.0 ) * 60.0

# Loop through results and count how many stars are in each bin
wmem = [x for x in range(len(mpri)) if mpri[x] > 0]
ntot = len(wmem)
nbins = np.zeros(len(mbins)-1)
for i in range(1,len(mbins)):
	wbin = [x for x in range(len(mpri)) if mpri[x] >= mbins[i][0] and mpri[x] < mbins[i][1]]
	nbins[i-1] = len(wbin)

# Compute radial results
radial = np.zeros([2, 10, len(mbins)])
totstep, binstep = int(ntot / len(radial[0,:,0])), int(max(nbins) / len(radial[0,:,0]))

print "\nUsing step size of", totstep, "for total distribution."
print "Using step size of", binstep, "for mass bin distribution."

for m in range(len(mbins)):
	# Detect stars in this mass range
	wbin = [x for x in range(len(mpri)) if mpri[x] > 0 and mpri[x] >= mbins[m][0] and mpri[x] < mbins[m][1]]
	if len(wbin) == 0: continue
	# Sort by radius
	rsrt = np.argsort([rad[x] for x in wbin])
	# Determine bin size
	if m == 0: stepsize = totstep
	else: stepsize = binstep
	# Loop through steps and compute statistics
	for s in range(len(radial[0,:,0])):
		if stepsize*(s+1) >= len(rsrt): continue
		wstep = rsrt[stepsize*s:stepsize*(s+1)]
		rstep = [rad[wbin[x]] for x in wstep]
		bstep = [x for x in wstep if msec[wbin[x]] > 0]
		radial[0,s,m] = np.mean(rstep)
		radial[1,s,m] = float(len(bstep)) / float(stepsize)
		
# Print results to file
of = open(options['data']+'--radial.txt', 'w')
for l in range(len(radial[0,:,0])):
	outstr = ''
	for m in range(len(mbins)): outstr += "%7.4f %5.3f  " % (radial[0,l,m] * 10.0**((options['m-M']+5.0)/5.0) * np.tan(1.0/60.0*np.pi/360.0), radial[1,l,m])
	print >>of, outstr
of.close()
		
# Plot Results
#pltcol = ['k-', 'b-', 'g-', 'c-', 'r-', 'y-']
#for m in range(len(mbins)):
#	plt.plot([x for x in radial[0,:,m] if x > 0], [x for x in radial[1,:,m] if x > 0], pltcol[m])
#plt.show()

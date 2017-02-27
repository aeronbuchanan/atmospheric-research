#!/usr/bin/python3
#
# Processing convective cloud data to compare SCAI and COP
#
# Data should be: 
# * one line per scene,
# * with comma-separated (x,y) position and radius for each convective cluster
# * convective clusters separated by semi-colons
# x0, y0, r0; x1, y1, r1; ...

import math
import argparse
import csv

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--inputdata", help="data csv file", default="data/data.csv")
parser.add_argument("-o", "--outputdir", help="directory to which to write output", default="results/")

args = parser.parse_args()
if args.outputdir[-1:] != '/': args.outputdir += '/'
print(args) # DEBUG

def dist(p1, p2):
  dx = p1[0] - p2[0]
  dy = p1[1] - p2[0]
  d = math.sqrt(dx * dx + dy * dy)
  return d

SQRTPI = math.sqrt(math.pi)
def dip(d, r1, r2):
  return (r1 + r2) / (d * SQRTPI)

Nmax = 200 # for SCAI calculation
L = 1000 # for SCAI calculation

odata = {}
sums = {}
counts = []
stats = ['arithmetic', 'geometric', 'COP', 'SCAI']
for s in stats:
  odata[s] = [False]
  sums[s] = [0.0]

maxN = 0

# load data
ln = 0
with open(args.inputdata) as f:
  for line in f:
    ln += 1

    cdata = []
    line = _s.rstrip().lstrip()
    for cs in csv.reader([line], quotechar='"', delimiter=';', quoting=csv.QUOTE_ALL, skipinitialspace=True):
      for vs in csv.reader(cs, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
        cdata.append(list(map(float, vs)))

    for s in stats: odata[s].append([False])
    while len(sums['arithmetic']) < len(cdata):
      counts.append(0)
      for s in stats: 
        sums[s].append(0.0)

    a = 0.0
    g = 1.0
    w = 0.0
    n = 0
    SCAI = 0.0
    COP = 0.0

    for i in range(1, len(cdata)):
      for j in range(0, i):
        n += 1
        d = dist(cdata[i][:2], cdata[j][:2])
        dd = dip(d, cdata[j][2], cdata[i][2])

        a += d
        g *= d
        w += dd

      # aggregate distance metrics for first (i + 1) clusters
      aa = a / n
      gg = pow(g, 1.0 / n)

      odata['arithmetic'][ln].append(aa)
      odata['geometric'][ln].append(gg)

      SCAI = i / Nmax * gg / L * 1000
      COP = w / n

      odata['SCAI'][ln].append(SCAI)
      odata['COP'][ln].append(COP)

      # collect sums for summary stats
      counts[i] += 1
      sums['arithmetic'][i] += aa
      sums['geometric'][i] += gg
      sums['SCAI'][i] += SCAI
      sums['COP'][i] += COP

    # save final values
    odata['arithmetic'][ln][0] = a / n
    odata['geometic'][ln][0] = pow(g, 1.0 / n)
    odata['SCAI'][ln][0] = SCAI
    odata['COP'][ln][0] = COP
    sums['arithmetic'][0] += aa
    sums['geometric'][0] += gg
    sums['SCAI'][0] += SCAI
    sums['COP'][0] += COP
    counts[0] += 1

    # TODO: check no overflows

    maxN = max(maxN, len(cdata))

outputdatafiles = []
correlations = []
for i in range(0, maxN): 
  correlations.append({c1: 0, c2: 0})
  outputdatafiles[].append(open(args.outputdir + 'output' + i + '.csv', 'w'))

means = {}
for s in stats: 
  means[s] = []
  for i in range(0, len(counts)):
    means[s].append( sums[s][i] / counts[i] )

variances = {}
for s in stats:
  variances[s] = []
  for i in range(0, len(counts)):
    variances[s].append(0.0)

covariances = {
  'arithgeo': variances['arithmetic'].copy(),
  'scaicop': variances['scai'].copy()
}

for i in range(1, ln):
  for j in range(0, len(odata['arithmetic'][i])):
    a = odata['arithmetic'][i][j]
    g = odata['geometric'][i][j]
    scai = odata['SCAI'][i][j]
    cop = odata['COP'][i][j]

    outputdatafiles[j].write('%f,%f,%f,%f\n' % (a, g, scai, cop))

    # stats
    for s in stats:
      d = odata[s][i][j] - means[s][j]
      variances[s][j] += d * d

    covariances['arithgeo'][j] += (a - means['arithmetic'][i][j]) * (g - means['geometric'][i][j])
    covariances['scaicop'][j] += (scai - means['scai'][i][j]) * (cop - means['cop'][i][j])

# final stats calc
for i in range(0, len(counts))
  covariances['arithgeo'][i] /= math.sqrt(variances['arithmetic'][i]) * math.sqrt(variances['geometric'][i])
  covariances['scaicop'][i] /= math.sqrt(variances['scai'][i]) * math.sqrt(variances['cop'][i])

for s in stats:
  for i in range(0, len(counts)):
    variances[s][i] /= (counts[i] - 1.0)

f = open(args.outputdir + 'correlations.csv', 'w')
for s in stats:
  f.write(s + " variances:\n")
  f.write(str(variances[s]) + "\n")
f.write("arith-geo covariances:\n")
f.write(str(covariances['arithgeo']) + "\n")
f.write("scai-cop covariance:\n")
f.write(str(covariances['scaicop']) + "\n")


















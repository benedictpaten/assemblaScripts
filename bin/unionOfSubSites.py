
import sys

"""Computes the union of sites where there is a substitution on the output from the substitution stats file.
"""

h = {}
for fH in [ open(f, 'r') for f in sys.argv[1:]]:
    for line in fH.readlines():
        if "INDEL-SUBSTITUTION:" in line:
            sequence, site = tuple(line.split()[1:])
            if h.has_key((sequence, site)):
                h[(sequence, site)] += 1
            else:
                h[(sequence, site)] = 1
    fH.close()

numberOfSites = sum(h.values())
numberOfCommonSites = len(h)

print numberOfSites, numberOfCommonSites
i = list(set(h.values()))
i.sort()
k = h.values()
k.sort()
for j in i:
    print j, k.count(j)


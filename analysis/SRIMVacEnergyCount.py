#!/usr/bin/python

#
# Tool to parse a SRIM COLLISION.txt file and produce an output comparable
# to TrimVacEnergyCount (output.type = 'vaccount')
#

import fileinput
import math
import re

recoil = re.compile('^\xdb')

# read file header
header = [''] * 4
for i in range(4) :
  header[i] = fileinput.input()

# parse rest of the file
vac = []
for line in fileinput.input() :
  # detect recoils
  if recoil.match(line) :
    field = line.split()

    # vacancy
    if field[7] == '1' :
      x = int(float(field[4]))
      E = int(math.log10(float(field[3])))
      E =  max(0, E);

      if E not in vac :
        vac += [[0]] * (1 + E - len(vac))

      try:
        vac[E][x] += 1
      except IndexError:
        vac[E] += [0] * (1 + x - len(vac[E]))
        vac[E][x] = 1

# output histogram
for E in range(len(vac)) :
  for x in range(len(vac[E])) :
    print "%d %d %d" % (E, x, vac[E][x])
  print

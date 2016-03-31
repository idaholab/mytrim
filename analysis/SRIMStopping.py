#!/usr/bin/python

#
# Tool to parse a XXXX on YYYY.txt stopping power table and produce an output
# comparable to runstopping
#

import fileinput
import sys

start = '  --------------  ---------- ---------- ----------  ----------  ----------'
end = '-----------------------------------------------------------'

table = False

# parse file
vac = []
for line in fileinput.input() :
  line = line.rstrip('\r\n');

  # detect end of data table
  if table :
    if line == end :
      table = False
      continue


  # detect start of data table
  if not table :
    if line == start :
      table = True
    continue

  # parse data
  field = line.split()
  if field[1] == 'eV' :
    mul = 1.0;
  elif field[1] == 'keV' :
    mul = 1000.0;
  elif field[1] == 'MeV' :
    mul = 1000000.0;
  else :
    sys.exit("Encountered invalid energy unit %s." % field[1])

  print "%f %f" % (float(field[0]) * mul, float(field[2]))

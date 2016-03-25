#!/usr/bin/python

#
# Tool to parse a SRIM COLLISION.txt file
#

import fileinput
import re

#pr = re.compile('Prime Recoil')
pr = re.compile('Recoil')

for rawline in fileinput.input() :
  line = re.sub(r'[^\w+-.]+', ' ', rawline).strip(' \t\n\r')

  # detect primary recoils
  if pr.match(line) :
    print line

#  else :
#    print "'%s'" % line

#!/usr/bin/python

import sys

f = open(sys.argv[0],'r')
while (read_data = f.readline()):
    print read_data

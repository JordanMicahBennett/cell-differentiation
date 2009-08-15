#!/usr/bin/python

import sys

## Read rule file
if len(sys.argv) > 1:
    f = open(sys.argv[1],'r')
    filedata = f.readlines()

## Get symbols
data = []
for x in range(0,len(filedata)):
    temp = filedata[x].split()
    for y in range(0,len(temp)):
        data.append(temp[y])
## Sort symbols
symbols = list(set(data))
symbols_inv = dict()
for x in range(0,len(symbols)):
    symbols_inv[symbols[x]] = x

## Get rules
rules = list()
for x in range(0,len(filedata)):
    temp = filedata[x].split()
    rules.append(list())
    for y in temp: rules[x].append(symbols_inv[y])


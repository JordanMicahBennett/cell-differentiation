#!/usr/bin/python

import sys

if len(sys.argv) != 4:
    print "Usage:",sys.argv[0],"<num_generations> <rule_file> <init_file>"
    sys.exit()

number_of_generations = int(sys.argv[1])
rule_file = sys.argv[2]
init_file = sys.argv[3]

if number_of_generations < 1:
    print "ERROR! Number of generations must be >= 1"
    sys.exit()

## Read rule file
f = open(rule_file,'r')
filedata = f.readlines()
f.close()

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

if len(rules) == 0:
    print "ERROR! No rules found in file:",rules_file
    sys.exit()

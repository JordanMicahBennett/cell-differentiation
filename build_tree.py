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
for x in range(len(filedata)):
    temp = filedata[x].split()
    for y in range(len(temp)):
        data.append(temp[y])
## Sort symbols
symbols = list(set(data))
symbols_inv = dict()
for x in range(len(symbols)):
    symbols_inv[symbols[x]] = x

## Get rules
rules = list()
for x in range(len(filedata)):
    temp = filedata[x].split()
    rules.append(list())
    for y in temp: rules[x].append(symbols_inv[y])

if len(rules) == 0:
    print "ERROR! No rules found in file:",rules_file
    sys.exit()

## Read initial state file
f = open(init_file,'r')
filedata = f.readlines()
f.close()

## Get initial states
init_states = []
try:
    for x in range(len(filedata)):
        init_states.append(list())
        temp = filedata[x].split()
        for y in temp: init_states[x].append(symbols_inv[y])
except:
    print "ERROR! Invalid symbols in initial state file."
    sys.exit()

if len(init_states) == 0:
    print "ERROR! No initial states found in file:",init_file


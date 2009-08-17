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
rule_list = []
for x in range(len(filedata)):
    temp = filedata[x].split()
    rule_list.append(list())
    for y in temp: rule_list[x].append(symbols_inv[y])

if len(rule_list) == 0:
    print "ERROR! No rules found in file:",rules_file
    sys.exit()

## Put rules in dictionary
rules = dict()
for x in range(len(symbols)):
    rules[x] = list()
    for y in rule_list: 
        if y[0] == x: rules[x].append(y[1:len(y)])

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
    print "ERROR! Invalid symbol in initial state file:",y
    sys.exit()

if len(init_states) == 0:
    print "ERROR! No initial states found in file:",init_file

## Node data structure
class Node:
    state = None
    generation = None
    expand = None
    def __init__(self,state,generation,expand):
        self.state = state
        self.generation = generation
        self.expand = expand
        
## Function to perform expansion
def expand(stack,final,number_of_generations,rules):
    n = stack.pop()
    if len(n.expand) == 0:
        if n.generation = number_of_generations:
            final.append(n.state)
        else:
            n.expand = range(len(n.state))
            n.generation = n.generation + 1
            stack.append(n)
    else:
        symbol_number = n.expand.pop()
        symbol = n.state[symbol_number]
        try:
            for x in range(len(rules[symbol])):
                stack.append(Node(n.state[:symbol_number] + x + 
                                  n.state[symbol_number+1:],
                                  n.generation,
                                  n.expand))
        except: ## No rule found for symbol (assume same for next generation)
            stack.append(Node(n.state,n.generation,n.expand))
    del n

## Initialize stack and final list
stack = list()
final = list()
for x in init_states:
    stack.append(Node(x,1,range(len(x))))

## Perform expansion
while len(stack) > 0: expand(stack,final,number_of_generations,rules)

print final

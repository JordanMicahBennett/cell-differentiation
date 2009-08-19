#!/usr/bin/python

import sys

debug = sys.stdout
##debug = open("/dev/null",'w')

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
    temp = filedata[x].split(':')[0].split()
    for y in temp:
        for z in y.split('*'):
            if not z.isdigit(): data.append(z)
## Sort symbols
symbols = list(set(data))
symbols_inv = dict()
for x in range(len(symbols)):
    symbols_inv[symbols[x]] = x
del data

## Get rules
rules = []
for x in range(len(filedata)):
    temp = filedata[x].split(':')[0].split()
    rules.append([])
    for y in temp:
        subtemp = y.split('*')
        if subtemp[0].isdigit():
            rules[x] = rules[x] + [symbols_inv[subtemp[1]]]*int(subtemp[0])
        else:
            rules[x].append(symbols_inv[y])

if len(rules) == 0:
    print "ERROR! No rules found in file:",rules_file
    sys.exit()

## Put rules in dictionary
rules_inv = dict()
for x in range(len(symbols)):
    rules_inv[x] = []
    for y in range(len(rules)):
        if rules[y][0] == x: rules_inv[x].append(y)
    if len(rules_inv[x]) == 0:
        del rules_inv[x]

## Read initial state file
f = open(init_file,'r')
filedata = f.readlines()
f.close()

## Get initial states
init_states = []
try:
    for x in range(len(filedata)):
        init_states.append([])
        temp = filedata[x].split(':')[0].split()
        for y in temp:
            subtemp = y.split('*')
            if subtemp[0].isdigit(): 
                init_states[x] = init_states[x] + [symbols_inv[subtemp[1]]]*int(subtemp[0])
            else:
                init_states[x].append(symbols_inv[y])
except:
    print "ERROR! Invalid symbol in initial state file:",y
    sys.exit()
init_states.reverse()

if len(init_states) == 0:
    print "ERROR! No initial states found in file:",init_file

## Node data structure
class Node:
    state = None
    generation = None
    expand = None
    selected = None
    def __init__(self,state,generation,expand,selected):
        self.state = state
        self.generation = generation
        self.expand = expand
        self.selected = selected
    def write(self,stream=sys.stdout):
        print >> stream,"State:",self.state
        print >> stream,"Generation:",self.generation
        print >> stream,"Need Expansion:",self.expand
        print >> stream,"Rules Selected:",self.selected
    def copy(self):
        return Node(list(self.state),
                    self.generation,
                    list(self.expand),
                    list(self.selected))

## Function to print results
def symbolize(state,symbols):
    return None

## Function to perform expansion
def expand(stack,final_set,number_of_generations,rules,rules_inv):
    print >> debug,"Stack:"
    for s in stack:
        s.write(debug)
    print >> debug
    n = stack.pop()
    print >> debug,"Starting:",n.state
    if len(n.expand) == 0:
        n.generation += 1
        final_set[n.generation].append(n.copy())
        if n.generation < number_of_generations:
            n.expand = range(len(n.state))
            stack.append(n)
        print >> debug
    else:
        symbol_number = n.expand.pop()
        print >> debug,"Modifying:",symbol_number
        try:
            for x in rules_inv[n.state[symbol_number]]:
                print >> debug,"First:",n.state[:symbol_number]
                print >> debug,"Middle:",rules[x][1:]
                print >> debug,"Last:",n.state[symbol_number+1:]
                print >> debug,"Ending:",n.state[:symbol_number] + rules[x][1:] + n.state[symbol_number+1:]
                new_node = n.copy()
                new_node.state = new_node.state[:symbol_number] + rules[x][1:] + new_node.state[symbol_number+1:]
                new_node.selected[x] += 1
                stack.append(new_node)
        except:
            print >> debug,"WARNING! No rules apply..."
            stack.append(n)
        print >> debug
    del n

## Initialize stack and final list
stack = []
final_set = []
for x in range(number_of_generations+1):
    final_set.append([])
for x in init_states:
    stack.append(Node(x,0,range(len(x)),[0]*len(rules)))
    final_set[0].append(stack[len(stack)-1].copy())

## Perform expansion
while len(stack) > 0: expand(stack,final_set,number_of_generations,rules,rules_inv)

## Create probability tables...
for n in range(len(final_set)):
    final = final_set[n]

    ## Sort state symbols
    for x in range(len(final)):
        final[x].state.sort()

    ## Sort final states
    final_states = []
    final_states_indexed = []
    for x in range(len(final)):
        final_states.append(str(final[x].state))
    final_states = list(set(final_states))
    final_states_probabilities = []
    for x in range(len(final_states)):
        final_states_probabilities.append([])
        final_states_indexed.append([])
        for y in range(len(final)):
            if str(final[y].state) == final_states[x]: 
                final_states_probabilities[x].append(list(final[y].selected))
                final_states_indexed[x] = final[y].state

    final_states = []
    for x in range(len(final_states_indexed)):
        final_states.append([0]*len(symbols))
        for y in range(len(symbols)):
            final_states[x][y] += final_states_indexed[x].count(y)

    print >> debug,"Generation",n
    print >> debug,final_states
    print >> debug,final_states_probabilities
    for x in range(len(final)):
        final[x].write(debug)
    print >> debug

    ## Output results
    filename = "generation_%03d.txt"%(n)
    f = open(filename,'w')

    for x in range(len(final_states)):
        for y in range(len(final_states[x])):
            if final_states[x][y]:
                if final_states[x][y] > 1:
                    print >> f,"%d*%s"%(final_states[x][y],symbols[y]),
                else:
                    print >> f,symbols[y],
        ## Print probabilities
        print >> f,":",

#         for y in range(len(rules)):
#             probabilities = dict()
#             for z in range(len(final_states_probabilities[x])):
#                 if final_states_probabilities[x][z][y] > 0:
#                     try:  
#                         probabilities[final_states_probabilities[x][z][y]] += 1
#                     except:
#                         probabilities[final_states_probabilities[x][z][y]] = 1    
#             for z in probabilities.items():
#                 if z[0] > 1:
#                     if z[1] > 1:
#                         print >> f,"%d*P%d^%d"%(z[1],y,z[0]),
#                     else:
#                         print >> f,"P%d^%d"%(y,z[0]),
#                 else:
#                     if z[1] > 1:
#                         print >> f,"%d*P%d"%(z[1],y),
#                     else:
#                         print >> f,"P%d"%(y),

        print >> f

    f.close()

## Extra testing stuff
print >> debug,"Symbols:"
print >> debug,symbols
print >> debug,"Inversion:"
print >> debug,symbols_inv
print >> debug,"Rules:"
print >> debug,rules
print >> debug,"Inversion:"
print >> debug,rules_inv

debug.close()

#!/usr/bin/python -u

from sympy import simplify
import sys
import gc
import time

##debug = sys.stdout
##debug = open("/dev/null",'w')

use_simplify = False

if len(sys.argv) != 4 and len(sys.argv) !=5 :
    print
    print "Usage:",sys.argv[0]," [-s] <num_generations> <rule_file> <init_file>"
    print "    -s : simplify the probabilities symbolically"
    print
    sys.exit()

if len(sys.argv) == 5:
    if sys.argv[1] == "-s":
        sys.argv.remove("-s")
        use_simplify = True
    else:
        print
        print "Usage:",sys.argv[0]," [-s] <num_generations> <rule_file> <init_file>"
        print "    -s : simplify the probabilities symbolically"
        print
        sys.exit()
        
number_of_generations = int(sys.argv[1])
rule_file = sys.argv[2]
init_file = sys.argv[3]

if number_of_generations < 1:
    print "ERROR! Number of generations must be >= 1"
    sys.exit()

## Minimal checking finished --- we are a go!

## Node data structure
class Node:
    state = None
    generation = None
    expand = None
    selected = None
    base_prob = None
    def __init__(self,state,generation,expand,selected,base_prob):
        self.state = state
        self.generation = generation
        self.expand = expand
        self.selected = selected
        self.base_prob = base_prob
    def write(self,stream=sys.stdout):
        print >> stream,"State:",self.state
        print >> stream,"Generation:",self.generation
        print >> stream,"Need Expansion:",self.expand
        print >> stream,"Rules Selected:",self.selected
        print >> stream,"Base Probability:",self.base_prob
    def copy(self):
        return Node(list(self.state),
                    self.generation,
                    list(self.expand),
                    list(self.selected),
                    self.base_prob)

## Function to perform expansion for one generation
def expand(stack,final,rules,rules_inv):
    n = stack.pop()
    if len(n.expand) == 0:
        n.generation += 1
        final.append(n.copy())
    else:
        symbol_number = n.expand.pop()
        try:
            for x in rules_inv[n.state[symbol_number]]:
                new_node = n.copy()
                new_node.state = new_node.state[:symbol_number] + rules[x][1:] + new_node.state[symbol_number+1:]
                new_node.selected[x] += 1
                stack.append(new_node)
        except:
            stack.append(n)
    del n

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
rule_probs = []
for x in range(len(filedata)):
    temp = filedata[x].split(':')
    rules.append([])
    try:
        rule_probs.append(temp[1].rstrip('\n'))
    except:
        rule_probs.append("P%d"%len(rules))
    temp = temp[0].split()
    for y in temp:
        subtemp = y.split('*')
        if subtemp[0].isdigit():
            rules[x] = rules[x] + [symbols_inv[subtemp[1]]]*int(subtemp[0])
        else:
            rules[x].append(symbols_inv[y])

if len(rules) == 0:
    print "ERROR! No rules found in file:",rules_file
    sys.exit()

## Are we using numbers or symbols for probabilities?
numeric = True
rule_probabilities = []
try:
    for x in range(len(rule_probs)):
        rule_probabilities.append(float(rule_probs[x]))
except:
    rule_probabilities = rule_probs
    numeric = False
del rule_probs

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
for x in range(len(filedata)):
    temp = filedata[x].split(':')
    state = [0] * len(symbols)
    for y in temp[0]:
        try:
            subtemp = y.split('*')
            if subtemp[0].isdigit(): 
                state[symbols_inv[subtemp[1]]] += int(subtemp[0])
            else:
                state[symbols_inv[y]] += 1
        except:
            print "ERROR! Invalid symbol in initial state file:",y
            sys.exit()
    
init_states.reverse()

if len(init_states) == 0:
    print "ERROR! No initial states found in file:",init_file

## Initialize stack and final list
stack = []
final = []
for x in init_states:
    final.append(Node(x,0,range(len(x)),[0]*len(rules)))

init_time = time.time()
## Perform expansion
for n in range(1,number_of_generations+1):
    print "Processing Generation",n
    print

    gen_start = time.time()

    for x in final:
        new_node = x.copy()
        x.expand = range(len(x.state))
        stack.append(x)
    final = []
    gc.collect()
    calls = 0
    while len(stack) > 0: 
        expand(stack,final,rules,rules_inv)
        calls += 1

    gen_end = time.time()
    print "Time elapsed:",(gen_end - gen_start)
    print "Expand function called",calls,"times."

    gen_start = time.time()

    ## Create probability tables and state representation
    final_states = dict()
    for x in range(len(final)):
        temp = [0] * len(symbols)
        for y in final[x].state:
            temp[y] += 1
        temp_str = ""
        for y in range(len(temp)):
            if temp[y] > 0:
                temp_str += " %d*%s"%(temp[y],symbols[y])
        temp_str = temp_str[1:]
        try :
            final_states[temp_str].append(final[x].selected)
        except:
            final_states[temp_str] = [final[x].selected]

    ## Recast probabilities
    if not numeric or use_simplify:
        for x in final_states.items():
            prob_strings = dict()
            for y in range(len(x[1])):
                prob_string = ""
                for z in range(len(x[1][y])):
                    if x[1][y][z] > 0:
                        if x[1][y][z] == 1:
                            prob_string += " * %s"%(rule_probabilities[z])
                        else:
                            prob_string += " * %s**%d"%(rule_probabilities[z],x[1][y][z])
                try:
                    prob_strings[prob_string] += 1
                except:
                    prob_strings[prob_string] = 1
            prob_string = "0.0 "
            for y in prob_strings.items():
                prob_string += "+ (%d %s) "%(y[1],y[0])
            if (use_simplify):
                prob_string = str(simplify(prob_string))
            final_states[x[0]] = prob_string
    else:
        for x in final_states.items():
            prob = 0
            for y in range(len(x[1])):
                prob_sub = 1
                for z in range(len(x[1][y])):
                    if x[1][y][z] > 0:
                        prob_sub *= rule_probabilities[z]**x[1][y][z]
                prob += prob_sub
            final_states[x[0]] = str(prob)

    ## Output results
    filename = "generation_%03d.txt"%(n)
    f = open(filename,'w')

    for x in final_states.items():
        print >> f,x[0],":",x[1]

    f.close()

    del final_states

    gen_end = time.time()
    print "Post-processing time:",(gen_end - gen_start)
    print

end_time = time.time()
print "Total elapsed time:",(end_time - init_time)

#debug.close()

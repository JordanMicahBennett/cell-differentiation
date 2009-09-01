#!/usr/bin/python -u

from sympy import simplify
from collections import deque
import os
import sys
import getopt
import gc
import time
import shelve

def usage():
    print
    print "Usage:",sys.argv[0]," [-s] <num_generations> <rule_file> <init_file>"
    print "    -s | --symbolic : simplify the probabilities symbolically"
    print
    sys.exit(-1)

use_simplify = False
use_threads = False
number_of_threads = 1

## Parse options
try:
    opts, args = getopt.getopt(sys.argv[1:],"s",["symbolic"])
except getopt.GetoptError:
    usage()

for opt,arg in opts:
    if opt in ("-s","--symbolic"):
        use_simplify = True
    elif opt in ("-t","--threads"):
        use_threads = True
        try:
            number_of_threads = int(arg)
        except:
            sys.stderr.write("\nERROR! Provided threads not an integer: %s\n\n"%(arg))
            sys.exit(-1)

if use_threads and number_of_threads < 1:
    sys.stderr.write("\nERROR! Number of threads must be greater than zero.\n\n")
    sys.exit(-1)
    
## Check for correct number of arguments
if len(args) != 3:
    usage()
        
number_of_generations = 0
rule_file = args[1]
init_file = args[2]

try :
    number_of_generations = int(args[0])
except:
    sys.stderr.write("\nERROR! Provided generations not an integer: %s\n\n"%(sys.argv[1]))
    usage()

if number_of_generations <= 0:
    sys.stderr.write("\nERROR! Number of generations must be greater than zero.\n\n")
    sys.exit(-1)

if not os.path.isfile(rule_file):
    sys.stderr.write("\nERROR! Rules file does not exist: %s\n\n"%(rule_file))
    sys.exit(-1)

if not os.path.isfile(init_file):
    sys.stderr.write("\nERROR! Initial state file does not exist: %s\n\n"%(init_file))
    sys.exit(-1)

## Minimal checking finished --- we are a go!

## Symbol table
## This class holds all information on the symbols, rules, and probabilities
## used during evaluation.
class SymbolTable:
    symbols = []
    symbols_inv = dict()
    rules = []
    rules_inv = dict()
    rules_probabilities = []
    use_numeric = True
    generic_symbol_number = 1
    def read_rules(self,rule_file):
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
        self.symbols = list(set(data))

        if len(self.symbols) == 0:
            sys.stderr.write("\nERROR! No symbols found in file: %s\n\n"%rule_file)
            return False

        ## Symbol inversion table
        for x in range(len(self.symbols)):
            self.symbols_inv[self.symbols[x]] = x
        del data

        ## Get rules
        rule_probs = []
        for x in range(len(filedata)):
            temp = filedata[x].split(':')
            self.rules.append([0] * len(self.symbols))
            try:
                prob = temp[1].replace(' ','').rstrip('\n')
                if prob:
                    rule_probs.append(prob)
                else:
                    rule_probs.append("P%d"%len(self.rules))
            except:
                rule_probs.append("P%d"%len(self.rules))
            temp = temp[0].split()
            first = -1
            for y in temp:
                subtemp = y.split('*')
                if subtemp[0].isdigit():
                    self.rules[x][self.symbols_inv[subtemp[1]]] += int(subtemp[0])
                    if first < 0:
                        first = self.symbols_inv[subtemp[1]]
                else:
                    self.rules[x][self.symbols_inv[subtemp[0]]] += 1
                    if first < 0:
                        first = self.symbols_inv[subtemp[0]]
            self.rules[x][first] -= 1
            try:
                self.rules_inv[first].append(x)
            except:
                self.rules_inv[first] = [x]
                    
        ## Are we using numbers or self.symbols for probabilities?
        self.use_numeric = True
        try:
            for x in range(len(rule_probs)):
                self.rules_probabilities.append(float(rule_probs[x]))
        except:
            self.rules_probabilities = rule_probs
            self.use_numeric = False
            del rule_probs

        ## Success!!
        return True

## Node data structure
## A single node of the tree is stored in this structure.
## It contains the counts of each symbol in the state it represents,
## the counts of each symbol that it needs to expand, the counts of each rule
## that was selected during expansion, the probability drawn from the initial
## state (or from the previous generation for continued processing,) and
## a reference to the symbol table used for evaluation.
class Node:
    state = None
    expandable = None
    selected = None
    base_prob = 1
    symbol_table = None
    def __init__(self,symbol_table):
        self.symbol_table = symbol_table
    def copy(self):
        new_node = Node(self.symbol_table)
        new_node.state = list(self.state)
        new_node.expandable = list(self.expandable)
        new_node.selected = list(self.selected)
        new_node.base_prob = self.base_prob
        return new_node
    def write(self,f):
        print >> f,"Object:",self
        print >> f,"State:",self.state
        print >> f,"Expandable:",self.expandable
        print >> f,"Selected:",self.selected
        print >> f,"Base Probability:",self.base_prob
        print >> f,"Symbol Table:",self.symbol_table
        print >> f,"As String:",self.tostring()
        print >> f
    def tostring(self):
        result = ""
        for x in range(len(self.state)):
            if self.state[x] > 0:
                if self.state[x] > 1:
                    result += "%d*%s "%(self.state[x],self.symbol_table.symbols[x])
                else:
                    result += "%s "%(self.symbol_table.symbols[x])
        result += ": "
        if (self.symbol_table.use_numeric):
            temp = self.base_prob
            for x in range(len(self.selected)):
                if self.selected[x] > 0:
                    temp *= pow(self.symbol_table.rules_probabilities[x],self.selected[x])
            result += str(temp)
        else:
            result += "%s"%(self.base_prob)
            for x in range(len(self.selected)):
                if self.selected[x] > 0:
                    if self.selected[x] > 1:
                        result += "*%s**%d"%(self.symbol_table.rules_probabilities[x],self.selected[x])
                    else:
                        result += "*%s"%(self.symbol_table.rules_probabilities[x])
        return result
    def fromstring(self,input=None):
        if input == None:
            return self
        temp = input.rstrip('\n').split(':')
        self.state = [0] * len(self.symbol_table.symbols)
        occupied = False
        for y in temp[0].split():
            try:
                subtemp = y.split('*')
                if subtemp[0].isdigit(): 
                    self.state[self.symbol_table.symbols_inv[subtemp[1]]] += int(subtemp[0])
                else:
                    self.state[self.symbol_table.symbols_inv[y]] += 1
                occupied = True
            except:
                sys.stderr.write("\nERROR! Invalid symbol in initial state spec: %s\n\n"%y)
                self.state = None
                return None
        if not occupied:
            self.state = None
            return None
        try:
            prob = temp[1].rstrip('\n').replace(' ','')
        except:
            prob = 1
        try:
            self.base_prob = float(prob)
        except:
            self.base_prob = str(prob)
            self.symbol_table.use_numeric = False
        self.expandable = list(self.state)
        self.selected = [0] * len(self.symbol_table.rules)
        return self                
    def expand(self,stack,shelf):
        expand = -1
        for x in range(len(self.expandable)):
            if self.expandable[x] > 0:
                expand = x
                break
        if expand >= 0:
            self.expandable[expand] -= 1
            self.state[expand] -= 1
            try:
                for x in self.symbol_table.rules_inv[expand]: 
                    n = self.copy()
                    for y in range(len(n.state)):
                        n.state[y] += self.symbol_table.rules[x][y]
                    n.selected[x] += 1
                    stack.append(n)
            except:
                self.state[expand] += 1
                self.expandable[expand] = 0
                stack.append(self)
        else:
            representation = self.tostring().split(" : ")
            if symbol_table.use_numeric:
                try:
                    prob = shelf[representation[0]]
                    prob += float(representation[1])
                    shelf[representation[0]] = prob
                except:
                    shelf[representation[0]] = float(representation[1])
            else:
                try:
                    prob_dict = shelf[representation[0]]
                    try:
                        prob_dict[representation[1]] += 1
                    except:
                        prob_dict[representation[1]] = 1
                except:
                    prob_dict = { representation[1] : 1 }
                shelf[representation[0]] = prob_dict
        return None

## Puts one state on the stack via a nodes read from a state file
def populate_stack(stack,symbol_table,infile):
    x = infile.readline()
    new_node = Node(symbol_table)
    new_node = new_node.fromstring(x)
    if new_node:
        stack.append(new_node)
        return True
    return False

## Puts all states into a list via a nodes read from a state file
def populate_list(symbol_table,infile):
    nodes = []
    in_f = open(infile,'r')
    for x in in_f:
        new_node = Node(symbol_table)
        new_node = new_node.fromstring(x)
        if new_node:
            nodes.append(new_node)
    in_f.close()
    return nodes

## Outputs the probabilities of each symbol ocurring
def symbol_count(states,symbol_table,outfile):
    if len(states) == 0:
        sys.stderr.write("\n\nERROR! No states generated from last generation!\n\n")
        sys.exit()
    out_f = open(outfile,'w')
    max_count = 0
    current_count = 0        
    ## Header info
    for x in symbol_table.symbols:
        print >> out_f,"\"%s\""%(x),
    print >> out_f
    ## Hit it!
    while current_count <= max_count:
        print >> out_f,current_count,
        for x in range(len(symbol_table.symbols)):
            if symbol_table.use_numeric:
                prob = 0
                for y in states:
                    if y.state[x] > max_count:
                        max_count = y.state[x]
                    if y.state[x] == current_count:
                        prob += y.base_prob
                print >> out_f,prob,
            else:
                prob = "0.0"
                for y in states:
                    if y.state[x] > max_count:
                        max_count = y.state[x]
                    if y.state[x] == current_count:
                        prob += " + %s"%(y.base_prob)
                if (use_simplify):
                    prob = str(simplify(prob)).replace(' ','')
                print >> out_f,"\"%s\""%(prob.replace("**","^")),
        print >> out_f
        current_count += 1
    out_f.close()
    return

## Reads in each item in a dictionay state:dict and combines the probabilities
## (arithmetic addition) of all states with the same number and type of symbols.
## In particular, this function works only on numeric probabilities.
def simplify_states_numeric(shelf,outfile):
    for state,probability in shelf.items():
        print >> outfile,state,":",probability
    
## Same as above function, but it works solely on symbolic probabilities. The
## flag "use_simplify" indicates whether to use Sympy to simplify the expressions
## on each generation (much slower, but probabilities are more compact and final
## evaluation would be faster.)
def simplify_states_symbolic(shelf,outfile,use_simplify):
    for state,value_dict in shelf.items():
        prob_dict = value_dict
        probability = "(0.0"
        for key_prob,value_count in prob_dict.items():
            if value_count > 1:
                probability += "+%d*(%s)"%(value_count,key_prob)
            else:
                probability += "+(%s)"%(key_prob)
        probability += ")"
        if (use_simplify):
            probability = "(" + str(simplify(probability)) + ")"
        print >> outfile,state,":",probability

## Functions and data structures defined... let us begin.

## Create symbol table and stack
symbol_table = SymbolTable()
stack = deque()

## Read rules into symbol_table
if not symbol_table.read_rules(rule_file):
    sys.exit(-1)

## Extra line
print

init_time = time.time()
## Perform expansion
for n in range(number_of_generations):
    print "Processing Generation",n+1
    print

    ## Drop garbage before this generation
    gc.collect()

    ## Read previous state file (or initial state file provided) and
    ## fill the stack with the states one at a time (expanding into
    ## temporary file.)
    gen_start = time.time()
    print "Expanding tree...",

    gen_shelf = shelve.open(".generation_%03d.%d.dat"%(n+1,os.getpid()))
    state_f = open(init_file,'r')
    calls = 0
    while populate_stack(stack,symbol_table,state_f):
        while len(stack) > 0: 
            stack.pop().expand(stack,gen_shelf)
            calls += 1

    if calls == 0:
        sys.stderr.write("\nERROR! Could not get states from file: %s\n\n"%(init_file))
        gen_shelf.close()
        os.remove(".generation.%03d.dat"%(n))
        sys.exit(-1)
        
    gen_end = time.time()
    print "done."
    print "Time elapsed:",(gen_end - gen_start)
    print "Expand function called",calls,"times."

    print "Post-processing state information...",
    gen_start = time.time()

    ## Catalog the current generation
    gen_file = open("generation_%03d.txt"%(n+1),'w')

    if symbol_table.use_numeric:
        simplify_states_numeric(gen_shelf,gen_file)
    else:
        simplify_states_symbolic(gen_shelf,gen_file,use_simplify)

    gen_file.close()
    gen_shelf.close()
    os.remove(".generation_%03d.%d.dat"%(n+1,os.getpid()))

    init_file = "generation_%03d.txt"%(n+1)

    ## Create summary table
    symbol_count(populate_list(symbol_table,init_file),symbol_table,"generation_%03d_summary.txt"%(n+1))

    gen_end = time.time()
    print "done."
    print "Time elapsed:",(gen_end - gen_start)
    print

end_time = time.time()
print "Total elapsed time:",(end_time - init_time)

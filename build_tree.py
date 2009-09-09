#!/usr/bin/python -u

from sympy import simplify
from collections import deque
import os
import sys
import getopt
import gc
import time
import shelve
import cPickle

def usage():
    print
    print 'Usage:',sys.argv[0],' [-s] [-e <epsilon>] <num_generations> <rule_file> <init_file>'
    print '    -s | --simplify : simplify the probabilities symbolically'
    print '    -e | --epsilon= : provide numerical cutoff for probabilities'
    print
    sys.exit(-1)

use_simplify = False
use_threads = False
number_of_threads = 1

## Calculate default epsilon - to machine precision
epsilon = 1
while epsilon / 2.0  + 1.0 > 1.0:
    epsilon = epsilon / 2.0

## Parse options
try:
    opts, args = getopt.getopt(sys.argv[1:],'se:',['simplify','epsilon='])
except getopt.GetoptError:
    usage()

for opt,arg in opts:
    if opt in ('-s','--simplify'):
        use_simplify = True
    elif opt in ('-t','--threads'):
        use_threads = True
        try:
            number_of_threads = int(arg)
        except:
            sys.stderr.write('\nERROR! Provided threads not an integer: %s\n\n'%(arg))
            sys.exit(-1)
    elif opt in ('-e','--epsilon'):
        try:
            new_epsilon = float(arg)
        except:
            sys.stderr.write('\nERROR! Provided epsilon not an real value: %s\n\n'%(arg))
            sys.exit(-1)
        if new_epsilon <= 0.0:
            sys.stderr.write('\nERROR! Epsilon must be a positive value.\n\n')
            sys.exit(-1)
            
        elif new_epsilon < epsilon:
            sys.stderr.write('\nWARNING! Provided epsilon is less than machine precision: %s\n\n'%(arg))
        elif new_epsilon >= 1.0:
            sys.stderr.write('\nERROR! Epsilon must be less than one.\n\n')
            sys.exit(-1)
        epsilon = new_epsilon
            
if use_threads and number_of_threads < 1:
    sys.stderr.write('\nERROR! Number of threads must be greater than zero.\n\n')
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
    sys.stderr.write('\nERROR! Provided generations not an integer: %s\n\n'%(sys.argv[1]))
    usage()

if number_of_generations <= 0:
    sys.stderr.write('\nERROR! Number of generations must be greater than zero.\n\n')
    sys.exit(-1)

if not os.path.isfile(rule_file):
    sys.stderr.write('\nERROR! Rules file does not exist: %s\n\n'%(rule_file))
    sys.exit(-1)

if not os.path.isfile(init_file):
    sys.stderr.write('\nERROR! Initial state file does not exist: %s\n\n'%(init_file))
    sys.exit(-1)

## Minimal checking finished --- we are a go!

def dump(object):
    return (cPickle.dumps(object,cPickle.HIGHEST_PROTOCOL))

def load(object):
    return (cPickle.loads(object))

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
            sys.stderr.write('\nERROR! No symbols found in file: %s\n\n'%rule_file)
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
                    rule_probs.append('('+prob+')')
                else:
                    rule_probs.append('P%d'%len(self.rules))
            except:
                rule_probs.append('P%d'%len(self.rules))
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
        for x in range(len(rule_probs)):
            temp = str(simplify(rule_probs[x]))
            try:
                self.rules_probabilities.append(float(temp))
            except:
                self.rules_probabilities.append('('+temp+')')
                self.use_numeric = False
        del rule_probs

        ## Success!!
        return True
    def state_to_string(self,state):
        result = ''
        for x in range(len(state)):
            if state[x] > 0:
                if state[x] > 1:
                    result += ' %d*%s'%(state[x],self.symbols[x])
                else:
                    result += ' %s'%(self.symbols[x])
        result = result[1:]
        return result
    def probability_to_string(self,prob,count=1):
        result = str(count)
        for x in range(len(prob)):
            if prob[x] > 0:
                if prob[x] > 1:
                    result += '*%s**%d'%(self.rules_probabilities[x],prob[x])
                else:
                    result += '*%s'%(self.rules_probabilities[x])
        return str(result)
    def probability_dict_to_string(self,prob_dict):
        result = '0.0'
        for prob,count in prob_dict.iteritems():
            prob = load(prob)
            if count > 0:
                result += '+' + self.probability_to_string(prob,count)
        return result
    def parse_state(self,input):
        temp = input.rstrip('\n').split(':')
        state = [0] * len(symbol_table.symbols)
        occupied = False
        for y in temp[0].split():
            try:
                subtemp = y.split('*')
                if subtemp[0].isdigit(): 
                    state[symbol_table.symbols_inv[subtemp[1]]] += int(subtemp[0])
                else:
                    state[symbol_table.symbols_inv[y]] += 1
                occupied = True
            except:
                sys.stderr.write('\nERROR! Invalid symbol in state spec: %s\n\n'%y)
                return None
        if not occupied:
            return None
        else:
            return (state)

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
    symbol_table = None
    def __init__(self,state,symbol_table):
        self.state = list(state)
        self.expandable = list(state)
        self.selected = [0] * len(symbol_table.rules)
        self.symbol_table = symbol_table
    def copy(self):
        new_node = Node(self.state,self.symbol_table)
        new_node.expandable = list(self.expandable)
        new_node.selected = list(self.selected)
        return new_node
    def write(self,f):
        print >> f,'Object:',self
        print >> f,'State:',self.state
        print >> f,'Expandable:',self.expandable
        print >> f,'Selected:',self.selected
        print >> f,'Symbol Table:',self.symbol_table
        print >> f,'As String:',self.tostring()
        print >> f
    def to_string(self):
        return (self.symbol_table.state_to_string(self.state) + ' : ' + 
                self.symbol_table.probability_to_string(self.selected),1)
    def expand(self,stack,shelf):
        expand = -1
        for x in range(len(self.expandable)):
            if self.expandable[x] > 0:
                expand = x
                break
        if expand >= 0:
            print
            print 'Intermediate',self.to_string()
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
                print
                print 'No Match!'
                self.state[expand] += 1
                self.expandable[expand] = 0
                stack.append(self)
        else:
            print
            print 'Final',self.to_string()
            try:
                prob_dict = load(shelf[dump(self.state)])
                try:
                    prob_dict[dump(self.selected)] += 1
                except:
                    prob_dict[dump(self.selected)] = 1
            except:
                prob_dict = { dump(self.selected) : 1 }
            shelf[dump(self.state)] = dump(prob_dict)
        return None

## Pulls some initial states from a file and populate the shelf
def populate_shelf(shelf,filename,symbol_table):
    f = open(filename,'r')
    for line in f:
        state = symbol_table.parse_state(line)
        if state:
            shelf[dump(state)] = dump({ dump([0] * len(symbol_table.rules)) : 1 })
    f.close()
    return

## Sums two probability vectors
def sum_prob(prob_one,prob_two):
    result = list(prob_one)
    for x in range(len(prob_two)):
        result[x] += prob_two[x]
    return result

## Multiplies two probability dictionaries
def multiply(base_prob_dict,prob_dict,result_dict):
    for new_prob,new_count in prob_dict.iteritems():
        new_prob = load(new_prob)
        for base_prob,base_count in base_prob_dict.iteritems():
            base_prob = load(base_prob)
            try:
                result_dict[dump(sum_prob(base_prob,new_prob))] += base_count * new_count
            except:
                result_dict[dump(sum_prob(base_prob,new_prob))] = base_count * new_count                

## Print states to a file
def print_states(shelf,symbol_table,filename):
    f = open(filename,'w')
    for state,prob_dict in shelf.iteritems():
        print >> f,symbol_table.state_to_string(load(state)),':',symbol_table.probability_dict_to_string(load(prob_dict))
    f.close()

## Functions and data structures defined... let us begin.

## Create symbol table and stack
symbol_table = SymbolTable()
stack = deque()
last_gen = dict()

## Read rules into symbol_table
if not symbol_table.read_rules(rule_file):
    sys.exit(-1)

## Read initial states
populate_shelf(last_gen,init_file,symbol_table)

## Extra line
print

init_time = time.time()
## Perform expansion
for n in range(number_of_generations):
    print 'Processing Generation',n+1
    print

    ## Drop garbage before this generation
    gc.collect()

    gen_start = time.time()

    ## Read previous state file (or initial state file provided) and
    ## fill the stack with the states one at a time (expanding into
    ## temporary file.)
    event_start = time.time()
    print 'Expanding tree...'

    gen_shelf = dict()

    state_file_size = len(last_gen) + 1
    state_file_count = 1
    state_f = open(init_file,'r')
    calls = 0
    for state,base_prob_dict in last_gen.iteritems():

        print 'Progress: %4.1f %% \r'%(state_file_count * (100.0 / state_file_size)),
        sys.stdout.flush()

        state = load(state)
        stack.append(Node(state,symbol_table))
        
        state_shelf = dict()

        while len(stack) > 0:
            stack.pop().expand(stack,state_shelf)
            calls += 1
        
        ## Filter results back to generation results by multiplication
        for new_state,prob_dict in state_shelf.iteritems():
            print
            print load(new_state)
            try:
                result = load(gen_shelf[new_state])
            except:
                result = dict()
            multiply(load(base_prob_dict),load(prob_dict),result)
            gen_shelf[new_state] = dump(result)    
                
        state_file_count += 1

    print 'Progress: %4.1f %% \r'%(100.0),

    if calls == 0:
        sys.stderr.write('\nERROR! Could not get states from file: %s\n\n'%(init_file))
        sys.exit(-1)
        
    event_end = time.time()
    print
    print 'Time elapsed:',(event_end - event_start)
    print 'Expand function called',calls,'times.'

    print 'Post-processing state information...'
    event_start = time.time()

    ## Catalog the current generation
    print_states(gen_shelf,symbol_table,'generation_%03d.txt'%(n+1))
    
    last_gen = gen_shelf

    event_end = time.time()
    print
    print 'Time elapsed:',(event_end - event_start)
    print 'Processing summary information...'
    event_start = time.time()

    ## Create summary table
#    symbol_count(populate_list(symbol_table,init_file),symbol_table,'generation_%03d_summary.txt'%(n+1))

    event_end = time.time()
    gen_end = time.time()
    print
    print 'Time elapsed:',(event_end - event_start)
    print 'Time for this generation:',(gen_end - gen_start)
    print

end_time = time.time()
print 'Total elapsed time:',(end_time - init_time)

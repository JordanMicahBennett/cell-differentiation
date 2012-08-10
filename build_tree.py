#!/usr/bin/python -u

from sympy import simplify
from sympy import Symbol
from collections import deque
from collections import defaultdict
from string import join
import os
import sys
import getopt
import gc
import time
import shelve
import cPickle
import glob

## Pickling - got any vinegar?
def dump(object):
    return (cPickle.dumps(object,cPickle.HIGHEST_PROTOCOL))
def load(object):
    return (cPickle.loads(object))

## Special dictionary factory routine
def zerodict():
    return (defaultdict(int))

## Root exit routine
def rootexit():
    sys.exit(-1)

def usage():
    print
    print 'Usage:',sys.argv[0],' <num_generations> <rule_file> <init_file>'
#    print '    -e | --epsilon= : provide numerical cutoff for probabilities'
    print
    rootexit()

## Symbol table
## This class holds all information on the symbols, rules, and probabilities
## used during evaluation.
class SymbolTable:
    symbols = []
    symbols_inv = dict()
    rules = []
    rules_inv = []
    rules_probabilities = []
    use_numeric = True
    def pack(self):
        data = []
        data.append(self.symbols)
        data.append(self.symbols_inv)
        data.append(self.rules)
        data.append(self.rules_inv)
        data.append(self.rules_probabilities)
        data.append(self.use_numeric)
        return data
    def unpack(self,data):
        self.symbols = data[0]
        self.symbols_inv = data[1]
        self.rules = data[2]
        self.rules_inv = data[3]
        self.rules_probabilities = data[4]
        self.use_numeric = data[5]
        return self
    def write(self,f=sys.stdout):
        print >> f,'Object:',self
        print >> f,'Symbols:',self.symbols
        print >> f,'Symbols Inv:',self.symbols_inv
        print >> f,'Rules:',self.rules
        print >> f,'Rules Inv:',self.rules_inv
        print >> f,'Rule Probabilities:',self.rules_probabilities
        print >> f
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

        # Prep rules_inv
        self.rules_inv = []
        for x in range(len(self.symbols)):
            self.rules_inv.append([])

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
            self.rules_inv[first].append(x)
                    
        ## Are we using numbers or self.symbols for probabilities?
        self.use_numeric = True
        for x in range(len(rule_probs)):
            temp = str(simplify(rule_probs[x])).replace(' ','')
            try:
                self.rules_probabilities.append(float(temp))
            except:
                self.rules_probabilities.append('('+temp+')')
                self.use_numeric = False
        del rule_probs

        ## Success!!
        return True
    def state_to_string(self,state):
        temp = []
        for x in range(len(state)):
            if state[x] > 0:
                if state[x] > 1:
                    temp.append('%d*%s'%(state[x],self.symbols[x]))
                else:
                    temp.append(self.symbols[x])
        return join(temp)
    def probability_to_string(self,prob,count=1):
        if count == 1:
            temp = []
        else:
            temp = [str(count)]
        for x in range(len(prob)):
            if prob[x] > 0:
                if prob[x] > 1:
                    temp.append('%s^%d'%(self.rules_probabilities[x],prob[x]))
                else:
                    temp.append(self.rules_probabilities[x])
        return join(temp,'*')
    def probability_to_c_string(self,prob,count=1):
        if count == 1:
            temp = []
        else:
            temp = ['%d.0'%(count)]
        for x in range(len(prob)):
            if prob[x] > 0:
                if prob[x] > 1:
                    temp.append('pow(rules[%d],%d.0)'%(x,prob[x]))
                else:
                    temp.append('rules[%d]'%(x))
        return join(temp,'*')
    def probability_dict_to_string(self,prob_dict):
        temp = []
        for prob,count in prob_dict.iteritems():
            prob = load(prob)
            if count > 0:
                temp.append(self.probability_to_string(prob,count))
        if len(temp) > 0:
            return join(temp,'+')
        return '0.0'
    def probability_dict_to_c_string(self,prob_dict,probabilities):
        temp = []
        for prob,count in prob_dict.iteritems():
            if count > 0:
                temp.append('result += %d.0*terms[%d];'%(count,probabilities[prob]))
        return join(temp,'\n')
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
    visited = None
    symbol_table = None
    def __init__(self,state,symbol_table):
        self.state = [0] * len(state)
        self.expandable = list(state)
        self.selected = [0] * len(symbol_table.rules)
        self.visited = False
        self.symbol_table = symbol_table
    def copy(self):
        new_node = Node(self.expandable,self.symbol_table)
        new_node.state = list(self.state)
        new_node.selected = list(self.selected)
        new_node.visited = self.visited
        return new_node
    def write(self,f=sys.stdout):
        print >> f,'Object:',self
        print >> f,'State:',self.state
        print >> f,'Expandable:',self.expandable
        print >> f,'Selected:',self.selected
        print >> f,'Visited:',self.visited
        print >> f,'Symbol Table:',self.symbol_table
    def expand(self,stack,leaf_dict,expand_shelf):
        if not self.visited:
            if not expand_shelf.has_key(dump(self.expandable)):
                expand_max = max(self.expandable)
                expand = self.expandable.index(expand_max)
                if expand_max > 0:
                    ## Put yourself back on the stack, but now labeled as visited
                    base = self.copy()
                    self.visited = True
                    stack.append(self)
                    base.expandable[expand] -= 1
                    base.state[expand] -= 1
                    for selected in self.symbol_table.rules_inv[expand]: 
                        n = base.copy()
                        for x in range(len(n.state)):
                            n.state[x] += self.symbol_table.rules[selected][x]
                        n.selected[selected] += 1
                        stack.append(n)
                else:
                    ## You are a leaf. You need to be in the leaf dictionary.
                    leaf_dict[dump(self.state)][dump(self.selected)] += 1
        else:
            new_dict = defaultdict(zerodict)
            expand_max = max(self.expandable)
            expand = self.expandable.index(expand_max)
            if expand_max > 0:
                expandable = list(self.expandable)
                expandable[expand] -= 1
                if expand_shelf.has_key(dump(expandable)):
                    subtree_dict = load(expand_shelf[dump(expandable)])
                    for selected in self.symbol_table.rules_inv[expand]:
                        state_update = list(self.symbol_table.rules[selected])
                        state_update[expand] -= 1
                        for state,prob_dict in subtree_dict.iteritems():
                            new_state = load(state)
                            for x in range(len(new_state)):
                                new_state[x] += state_update[x]
                            for prob,count in prob_dict.iteritems():
                                new_prob = load(prob)
                                new_prob[selected] += 1
                                new_dict[dump(new_state)][dump(new_prob)] += count
                else:
                    for state,prob_dict in leaf_dict.iteritems():
                        new_state = load(state)
                        for x in range(len(new_state)):
                            new_state[x] -= self.state[x]
                        for prob,count in prob_dict.iteritems():
                            new_prob = load(prob)
                            for x in range(len(new_prob)):
                                new_prob[x] -= self.selected[x]
                            new_dict[dump(new_state)][dump(new_prob)] += count
                    leaf_dict = defaultdict(zerodict)
            expand_shelf[dump(self.expandable)] = dump(new_dict)                    
        return None

## Pulls some initial states from a file and populate the shelf
def populate_shelf(shelf,filename,symbol_table):
    f = open(filename,'r')
    for line in f:
        state = symbol_table.parse_state(line)
        if state:
            shelf[dump(state)] = { dump([0] * len(symbol_table.rules)) : 1 }
    f.close()
    return

##Multiplies two probability dictionaries
def multiply(base_prob_dict,prob_dict,result_dict):
    for new_prob,new_count in prob_dict.iteritems():
        new_prob = load(new_prob)
        for base_prob,base_count in base_prob_dict.iteritems():
            base_prob = load(base_prob)
            temp = dump(add_prob_prob(base_prob,new_prob))
            if result_dict.has_key(temp):
                result_dict[temp] += base_count * new_count
            else:
                result_dict[temp] = base_count * new_count

def add_prob_prob(dest_prob,src_prob):
    result = list(dest_prob)
    for x in range(len(src_prob)):
        result[x] += src_prob[x]
    return result

## Integrate two shelf/dictionaries
def add_shelf_shelf(dest_shelf,src_shelf):
    for state,src_prob_dict in src_shelf.iteritems():
        add_shelf_prob_dict(dest_shelf,state,src_prob_dict)

def add_shelf_prob_dict(dest_shelf,state,src_prob_dict):
    if dest_shelf.has_key(state):
        dest_prob_dict = dest_shelf[state]
    else:
        dest_prob_dict = dict()
    add_prob_dict_prob_dict(dest_prob_dict,src_prob_dict)
    dest_shelf[state] = dest_prob_dict

def add_prob_dict_prob_dict(dest_prob_dict,src_prob_dict):
    for src_prob,src_count in src_prob_dict.iteritems():
        add_prob_dict_prob(dest_prob_dict,src_prob,src_count)

def add_prob_dict_prob(dest_prob_dict,src_prob,src_count):
    if dest_prob_dict.has_key(src_prob):
        dest_count = dest_prob_dict[src_prob]
    else:
        dest_count = 0
    dest_prob_dict[src_prob] = dest_count + src_count

## Print states to a file
def print_states(shelf,symbol_table,filename):
    out_f = open(filename,'w')
    max_count = len(shelf)+1
    current_count = 0
    for state,prob_dict in shelf.iteritems():
        state = load(state)
        print "Progress: %4.1f %% \r"%(current_count * (100.0 / max_count)),
        sys.stdout.flush()
        print >> out_f,symbol_table.state_to_string(state),':',symbol_table.probability_dict_to_string(prob_dict)
        current_count += 1
    print "Progress: %4.1f %% \r"%(100.0),
    sys.stdout.flush()
    out_f.close()
    return

## Print summary table
def print_c_code(states,probabilities,size,symbol_table,shelf,filename):

    ## Expressions
    for state,x in states.iteritems():
        out_f = open(filename + '_expression_%d.c'%(x),'w')
        print >> out_f,'void expression_%d(double* terms, double table[][%d]) {'%(x,len(symbol_table.symbols))
        print >> out_f,'double result = 0.0;'
        print >> out_f,symbol_table.probability_dict_to_c_string(shelf[state],probabilities)
        state = load(state)
        for y in range(len(state)):
            print >> out_f,'table[%d][%d] += result;'%(state[y],y)
        print >> out_f,'return;}'
        out_f.close()

    ## Terms
    terms_count = 0
    out_f = open(filename + '_terms_%d.c'%(terms_count / 200),'w')
    print >> out_f,'#include <math.h>'
    print >> out_f,'void terms_%d(double* rules, double* probs) {'%(0)
    for prob,x in probabilities.iteritems():
        terms_count += 1
        print >> out_f,'probs[%d] = %s;'%(x,symbol_table.probability_to_c_string(load(prob)))        
        if terms_count % 200 == 0:
            print >> out_f,'return;}'
            out_f.close()
            out_f = open(filename + '_terms_%d.c'%(terms_count / 200),'w')
            print >> out_f,'#include <math.h>'
            print >> out_f,'void terms_%d(double* rules, double* probs) {'%(terms_count / 200)
    print >> out_f,'return;}'
    out_f.close()

    ## Find unique symbols in rule probabilities
    unique_symbols = set()
    for rule_prob in symbol_table.rules_probabilities:
        unique_symbols = unique_symbols.union(simplify(rule_prob).atoms(Symbol))
    unique_symbols = list(unique_symbols)

    ## Probabilities
    out_f = open(filename + '_probabilities.c','w')
    print >> out_f,'#include <stdio.h>'
    print >> out_f,'#include <stdlib.h>'
    print >> out_f,'#include <malloc.h>'
    print >> out_f
    print >> out_f,'int main(int argc, char* argv[]) {'
    print >> out_f,'if (argc != %d) {'%(len(unique_symbols)+1)
    print >> out_f,'printf(\"\\nUsage: %s',
    for symbol in unique_symbols: print >> out_f,'%s'%(symbol),
    print >> out_f,'\\n\\n\",argv[0]);'
    print >> out_f,'return -1; }'
    print >> out_f,'double total;'
    for x in range(len(unique_symbols)):
        print >> out_f,'double %s = atof(argv[%d]);'%(unique_symbols[x],x+1)
    print >> out_f,'double* rules = malloc(sizeof(double)*%d);'%(len(symbol_table.rules_probabilities))
    for x in range(len(symbol_table.rules_probabilities)):
        print >> out_f,'rules[%d] = %s;'%(x,symbol_table.rules_probabilities[x]);
    print >> out_f,'printf(\"Symbol\\tRule\\tProbability\\n\");'
    for x in range(len(symbol_table.symbols)):
        print >> out_f,'total = 0.0;'
        for y in range(len(symbol_table.rules_inv[x])):
            print >> out_f,'printf(\"%s\\t%s'%(symbol_table.symbols[x],
                                            symbol_table.state_to_string(symbol_table.rules[symbol_table.rules_inv[x][y]])),
            print >> out_f,'\\t%0.18G\\n\",',
            print >> out_f,'rules[%d]);'%(symbol_table.rules_inv[x][y])
            print >> out_f,'total += rules[%d];'%(symbol_table.rules_inv[x][y])
        ## print >> out_f,'if (total < 1.0) {'
        print >> out_f,'   printf(\"%s\\t%s'%(symbol_table.symbols[x],
                                              symbol_table.symbols[x]),
        print >> out_f,'\\t%0.18G\\n\", 1.0-total);'
        ## print >> out_f,'}'

    print >> out_f,'free(rules);'
    print >> out_f,'return 0;}'
    out_f.close()

    ## Main
    out_f = open(filename + '.c','w')
    print >> out_f,'#include <stdio.h>'
    print >> out_f,'#include <stdlib.h>'
    print >> out_f,'#include <malloc.h>'
    print >> out_f
    print >> out_f,'// Expressions that compute the probabilities'
    for count in range(len(states)):
        print >> out_f,'void expression_%d(double*,double[][%d]);'%(count,len(symbol_table.symbols))
    print >> out_f
    print >> out_f,'// Precomputed terms for the probabilities'
    for count in range((len(probabilities)/200)+1):
        print >> out_f,'void terms_%d(double*,double*);'%(count)

    print >> out_f,'int main(int argc, char* argv[]) {'
    print >> out_f,'if (argc != %d) {'%(len(unique_symbols)+1)
    print >> out_f,'printf(\"\\nUsage: %s',
    for symbol in unique_symbols: print >> out_f,'%s'%(symbol),
    print >> out_f,'\\n\\n\",argv[0]);'
    print >> out_f,'return -1; }'
    for x in range(len(unique_symbols)):
        print >> out_f,'double %s = atof(argv[%d]);'%(unique_symbols[x],x+1)
    print >> out_f,'int x,y;'
    print >> out_f,'double table[%d][%d];'%(size,len(symbol_table.symbols))
    print >> out_f,'for (x = 0; x < %d; x++) for (y = 0; y < %d; y++) table[x][y] = 0;'%(size,len(symbol_table.symbols))
    print >> out_f,'double* rules = malloc(sizeof(double)*%d);'%(len(symbol_table.rules_probabilities))
    print >> out_f,'double* probs = malloc(sizeof(double)*%d);'%(len(probabilities))

    ## Pre-calculate those values
    for x in range(len(symbol_table.rules_probabilities)):
        print >> out_f,'rules[%d] = %s;'%(x,symbol_table.rules_probabilities[x]);
    for x in range((len(probabilities)/200)+1):
        print >> out_f,'terms_%d(rules,probs);'%(x)
    for state,x in states.iteritems():
        print >> out_f,'expression_%d(probs,table);'%(x)

    for x in symbol_table.symbols:
        print >> out_f,'printf(\"%s \");'%(x)
    print >> out_f,'printf(\"\\n\");'
    for x in range(size):
#        print >> out_f,'printf(\"%d \");'%x
        for y in range(len(symbol_table.symbols)):
            print >> out_f,'printf(\"%0.18G \",',
            print >> out_f,'table[%d][%d]);'%(x,y)
        print >> out_f,'printf(\"\\n\");'
    print >> out_f,'free(rules);'
    print >> out_f,'free(probs);'
    print >> out_f,'return 0; }'
    out_f.close()

    ## Makefile
    out_f = open('Makefile','a')
    print >> out_f,filename,':','%s_probabilities %s.o'%(filename,filename),
    for x in range(len(states)):
        print >> out_f,'%s_expression_%d.o'%(filename,x),
    for x in range((len(probabilities)/200)+1):
        print >> out_f,'%s_terms_%d.o'%(filename,x),
    print >> out_f
    print >> out_f,'\tfind \. -name \'%s_expression_*.o\' | xargs ar cr %s_expressions.a'%(filename,filename)
    print >> out_f,'\tfind \. -name \'%s_terms_*.o\' | xargs ar cr %s_terms.a'%(filename,filename)
    print >> out_f,'\tcc -lm -O3 -o %s'%(filename),'%s.o'%(filename),'%s_expressions.a'%(filename),'%s_terms.a'%(filename)
    print >> out_f
    print >> out_f,'%s_probabilities'%(filename),':','%s_probabilities.c'%(filename)
    print >> out_f,'\tcc -O3 -o %s_probabilities'%(filename),'%s_probabilities.c'%(filename)
    print >> out_f
    print >> out_f,'%s.o'%(filename),':','%s.c'%(filename)
    print >> out_f,'\tcc -O3 -c -o %s.o'%(filename),'%s.c'%(filename)
    for x in range(len(states)):
        print >> out_f,'%s_expression_%d.o'%(filename,x),':','%s_expression_%d.c'%(filename,x)
        print >> out_f,'\tcc -O3 -c -o %s_expression_%d.o'%(filename,x),'%s_expression_%d.c'%(filename,x)
    for x in range((len(probabilities)/200)+1):
        print >> out_f,'%s_terms_%d.o'%(filename,x),':','%s_terms_%d.c'%(filename,x)
        print >> out_f,'\tcc -O3 -c -o %s_terms_%d.o'%(filename,x),'%s_terms_%d.c'%(filename,x)
    print >> out_f
    print >> out_f
    out_f.close()

def expand_state(state,base_prob_dict,gen_shelf,expand_shelf,symbol_table):
    state_node = Node(state,symbol_table)
    for x in range(len(symbol_table.symbols)):
        if len(symbol_table.rules_inv[x]) == 0:
            state_node.state[x] += state_node.expandable[x]
            state_node.expandable[x] = 0

    if not expand_shelf.has_key(dump(state_node.expandable)):
        stack = deque()
        leaf_dict = defaultdict(zerodict)

        stack.append(state_node)
        while len(stack) > 0:
            stack.pop().expand(stack,leaf_dict,expand_shelf)

    if not expand_shelf.has_key(dump(state_node.expandable)):
        subtree = defaultdict(zerodict)
        for new_state,prob_dict in leaf_dict.iteritems():
            new_state = load(new_state)
            for x in range(len(new_state)):
                new_state[x] -= state[x]
            subtree[dump(new_state)] = prob_dict
    else:
        subtree = load(expand_shelf[dump(state_node.expandable)])
    
    ## Filter results back to generation results by multiplication
    for state_offset,prob_dict in subtree.iteritems():
        new_state = load(state_offset)
        for x in range(len(new_state)):
            new_state[x] += state[x]
        if gen_shelf.has_key(dump(new_state)):
            result = gen_shelf[dump(new_state)]
        else:
            result = dict()
        multiply(base_prob_dict,prob_dict,result)
        gen_shelf[dump(new_state)] = result

## Functions and data structures defined... let us begin.

symbol_table = SymbolTable()

## Calculate default epsilon - to machine precision
epsilon = 1
while epsilon / 2.0  + 1.0 > 1.0:
    epsilon = epsilon / 2.0

## Parse options
try:
    opts, args = getopt.getopt(sys.argv[1:],'e:',['simplify','epsilon='])
except getopt.GetoptError:
    usage()

for opt,arg in opts:
    if opt in ('-e','--epsilon'):
        try:
            new_epsilon = float(arg)
        except:
            sys.stderr.write('\nERROR! Provided epsilon not an real value: %s\n\n'%(arg))
            rootexit()
        if new_epsilon <= 0.0:
            sys.stderr.write('\nERROR! Epsilon must be a positive value.\n\n')
            rootexit()
        elif new_epsilon < epsilon:
            sys.stderr.write('\nWARNING! Provided epsilon is less than machine precision: %s\n\n'%(arg))
        elif new_epsilon >= 1.0:
            sys.stderr.write('\nERROR! Epsilon must be less than one.\n\n')
            rootexit()
        epsilon = new_epsilon

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
    rootexit()

if number_of_generations <= 0:
    sys.stderr.write('\nERROR! Number of generations must be greater than zero.\n\n')
    rootexit()

if not os.path.isfile(rule_file):
    sys.stderr.write('\nERROR! Rules file does not exist: %s\n\n'%(rule_file))
    rootexit()

if not os.path.isfile(init_file):
    sys.stderr.write('\nERROR! Initial state file does not exist: %s\n\n'%(init_file))
    rootexit()

## Minimal checking finished --- we are a go!

## Create last generation shelf
last_gen = shelve.open('.generation_%03d.%d'%(0,os.getpid()))

## Create state transition storage shelf
expand_shelf = shelve.open('.expand.%d'%(os.getpid()))

## Read rules into symbol_table
if not symbol_table.read_rules(rule_file):
    rootexit()

## Read initial states
populate_shelf(last_gen,init_file,symbol_table)

## Start new Makefile
f = open('Makefile','w')
print >> f,'all :',
for x in range(1,number_of_generations+1):
    print >> f,'generation_%03d_summary'%(x),
print >> f
print >> f
print >> f,'clean :'
for x in range(1,number_of_generations+1):
    print >> f,'\tfind \. -name \'generation_%03d_summary_expression_*.o\' | xargs rm'%(x)
    print >> f,'\tfind \. -name \'generation_%03d_summary_terms_*.o\' | xargs rm'%(x)
    print >> f,'\trm',
    print >> f,'generation_%03d_summary.o generation_%03d_summary_expressions.a generation_%03d_summary_terms.a'%(x,x,x)
print >> f
print >> f
f.close()

if len(last_gen) <= 0:
    sys.stderr.write('\nERROR! Initial state file has no data: %s\n\n'%(init_file))
    rootexit()

## Extra line
print

init_time = time.time()
## Perform expansion
for n in range(1,number_of_generations+1):
    print 'Processing Generation',n
    print

    gen_shelf = shelve.open('.generation_%03d.%d'%(n,os.getpid()))

    ## Drop garbage before this generation
    gc.collect()

    gen_start = time.time()

    ## Read previous state file (or initial state file provided) and
    ## fill the stack with the states one at a time (expanding into
    ## temporary file.)
    event_start = time.time()

    ## Send out work
    gen_size = len(last_gen) + 1
    gen_count = 1
    for state,base_prob_dict in last_gen.iteritems():

        gc.collect()

        print 'Expansion Progress: %4.1f %% \r'%(gen_count * (100.0 / gen_size)),
        sys.stdout.flush()

        expand_state(load(state),base_prob_dict,gen_shelf,expand_shelf,symbol_table)

        gen_count += 1

    print 'Expansion Progress: %4.1f %% \r'%(100.0)

    event_end = time.time()
    print 'Time elapsed:',(event_end - event_start)

    print 'Post-processing state information...'
    event_start = time.time()

    ## Catalog the current generation
    print_states(gen_shelf,symbol_table,'generation_%03d.txt'%(n))

    event_end = time.time()
    print
    print 'Time elapsed:',(event_end - event_start)
    print 'Generating C code...'
    event_start = time.time()

    ## Create summary table
    max_count = 0
    probabilities = shelve.open('.probabilities_%03d.%d'%(n,os.getpid()))
    states = shelve.open('.states_%03d.%d'%(n,os.getpid()))
    temp_count = 0
    for state in gen_shelf.iterkeys():
        states[state] = temp_count
        temp_count += 1
    temp_count = 0
    for state in states.iterkeys():
        for prob in gen_shelf[state].iterkeys():
            if not probabilities.has_key(prob):
                probabilities[prob] = temp_count
                temp_count += 1
        state = load(state)
        if max(state) > max_count:
            max_count = max(state)
    max_count += 1

    print_c_code(states,
                 probabilities,
                 max_count,
                 symbol_table,
                 gen_shelf,
                 'generation_%03d_summary'%(n))

    event_end = time.time()
    gen_end = time.time()
    print 'Time elapsed:',(event_end - event_start)
    print 'Time for this generation:',(gen_end - gen_start)
    print

    ## Promote to next generation
    last_gen.close()
    probabilities.close()
    states.close()
    for filename in glob.glob('.generation_%03d.%d*'%(n-1,os.getpid())):
        os.remove(filename)
    for filename in glob.glob('.probabilities_%03d.%d*'%(n,os.getpid())):
        os.remove(filename)
    for filename in glob.glob('.states_%03d.%d*'%(n,os.getpid())):
        os.remove(filename)
    last_gen = gen_shelf

end_time = time.time()
print 'Total elapsed time:',(end_time - init_time)

## Finalize results
last_gen.close()
expand_shelf.close()
for filename in glob.glob('.generation_%03d.%d*'%(number_of_generations,os.getpid())):
    os.remove(filename)
for filename in glob.glob('.expand.%d*'%(os.getpid())):
    os.remove(filename)

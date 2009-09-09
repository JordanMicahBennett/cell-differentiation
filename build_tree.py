#!/usr/bin/python -u

use_mpi = True
try:
    from mpi4py import MPI
except:
    use_mpi = False

from sympy import simplify
from collections import deque
import os
import sys
import getopt
import gc
import time
import shelve
import cPickle

mpi_rank = 0
mpi_size = 1
comm = None
if use_mpi:
    comm = MPI.COMM_WORLD
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()

## Pickling - got any vinegar?
def dump(object):
    return (cPickle.dumps(object,cPickle.HIGHEST_PROTOCOL))
def load(object):
    return (cPickle.loads(object))

## Root exit routine
def rootexit():
    comm.bcast(None,root=0)
    sys.exit(-1)

def usage():
    print
    print 'Usage:',sys.argv[0],' [-s] [-e <epsilon>] <num_generations> <rule_file> <init_file>'
    print '    -s | --simplify : simplify the probabilities symbolically'
    print '    -e | --epsilon= : provide numerical cutoff for probabilities'
    print
    rootexit()


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
    out_f = open(filename,'w')
    max_count = len(shelf)+1
    current_count = 0
    if use_simplify:
        for state,prob_dict in shelf.iteritems():
            print "Progress: %4.1f %% \r"%(current_count * (100.0 / max_count)),
            sys.stdout.flush()
            print >> out_f,symbol_table.state_to_string(load(state)),':',str(simplify(symbol_table.probability_dict_to_string(load(prob_dict)))).replace(' ','').replace('**','^')
            current_count += 1
    else:
        for state,prob_dict in shelf.iteritems():
            print "Progress: %4.1f %% \r"%(current_count * (100.0 / max_count)),
            sys.stdout.flush()
            print >> out_f,symbol_table.state_to_string(load(state)),':',symbol_table.probability_dict_to_string(load(prob_dict)).replace('**','^')
            current_count += 1
    print "Progress: %4.1f %% \r"%(100.0),
    sys.stdout.flush()
    out_f.close()
    return

## Print state summary information
def print_summary(shelf,symbol_table,filename):
    out_f = open(filename,'w')
    max_count = 0
    current_count = 0
    ## Header
    for x in symbol_table.symbols:
        print >> out_f,"\"%s\""%(x),
    print >> out_f
    # Hit it!
    while current_count <= max_count:
        print "Progress: %4.1f %% \r"%(current_count * (100.0 / (max_count + 1))),
        sys.stdout.flush()
        print >> out_f,current_count,
        ## Data structures
        sym_dict = []
        for x in range(len(symbol_table.symbols)):
            sym_dict.append(dict())
        for state,prob_dict in shelf.iteritems():
            state = load(state)
            prob_dict = load(prob_dict)
            for x in range(len(symbol_table.symbols)):
                if state[x] > max_count:
                    max_count = state[x]
                if state[x] == current_count:
                    for prob,count in prob_dict.iteritems():
                        try:
                            sym_dict[x][prob] += count
                        except:
                            sym_dict[x][prob] = count
        if use_simplify:
            for x in range(len(symbol_table.symbols)):
                print >> out_f,"\"%s\""%(str(simplify(symbol_table.probability_dict_to_string(sym_dict[x]))).replace(' ','').replace('**','^')),
        else:
            for x in range(len(symbol_table.symbols)):
                print >> out_f,"\"%s\""%(symbol_table.probability_dict_to_string(sym_dict[x]).replace('**','^')),
        print >> out_f
        current_count += 1
    out_f.close()
    print "Progress: %4.1f %% \r"%(100.0),
    sys.stdout.flush()
    return

## Functions and data structures defined... let us begin.

symbol_table = None

if mpi_rank == 0:
    use_simplify = False

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
        elif opt in ('-e','--epsilon'):
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

    ## Create symbol table and stack
    symbol_table = SymbolTable()
    last_gen = dict()

    ## Read rules into symbol_table
    if not symbol_table.read_rules(rule_file):
        rootexit()

    ## Read initial states
    populate_shelf(last_gen,init_file,symbol_table)

    if len(last_gen) <= 0:
        sys.stderr.write('\nERROR! Initial state file has no data: %s\n\n'%(init_file))
        rootexit()

    ## Extra line
    print

    ## Good to go.. broadcast the symbol table!
    symbol_table = comm.bcast(symbol_table,root=0)

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

        ## Send out work
        gen_size = len(last_gen) + 1
        gen_count = 1
        procs = set(range(1,mpi_size))
        working = set()
        for state,base_prob_dict in last_gen.iteritems():

            print 'Sending Progress: %4.1f %% \r'%(gen_count * (100.0 / gen_size)),
            sys.stdout.flush()

            dest = comm.recv(source=MPI.ANY_SOURCE,tag=1)
            comm.send('EXPAND',dest=dest,tag=2)
            comm.send(state,dest=dest,tag=3)
            comm.send(base_prob_dict,dest=dest,tag=4)

            working = working.union([dest])
            gen_count += 1

        print 'Sending Progress: %4.1f %% \r'%(100.0)

        sleeping = set()
        gen_size = mpi_size + 1
        gen_count = 1
        ## Gather all results
        while len(working) > 0:

            print 'Gathering Progress: %4.1f %% \r'%(gen_count * (100.0 / gen_size)),

            while 1:
                dest = comm.recv(source=MPI.ANY_SOURCE,tag=1)
                if dest in working:
                    break
                comm.send('WAIT',dest=dest,tag=2)
                sleeping = sleeping.union([dest])
                
            comm.send('COMBINE',dest=dest,tag=2)
            prob_dict = comm.recv(source=dest,tag=3)
            working.remove(dest)

            ## Integrate
            
            gen_count += 1

        for x in sleeping:
            comm.send('WAKEUP',dest=x,tag=3)

        print 'Gathering Progress: %4.1f %% \r'%(100.0)

#             state_shelf = dict()

#             while len(stack) > 0:
#                 stack.pop().expand(stack,state_shelf)
#                 calls += 1

#             ## Filter results back to generation results by multiplication
#             for new_state,prob_dict in state_shelf.iteritems():
#                 try:
#                     result = load(gen_shelf[new_state])
#                 except:
#                     result = dict()
#                 multiply(load(base_prob_dict),load(prob_dict),result)
#                 gen_shelf[new_state] = dump(result)    

#             state_count += 1

        event_end = time.time()
        print
        print 'Time elapsed:',(event_end - event_start)

        print 'Post-processing state information...'
        event_start = time.time()

        ## Catalog the current generation
        print_states(gen_shelf,symbol_table,'generation_%03d.txt'%(n+1))

        event_end = time.time()
        print
        print 'Time elapsed:',(event_end - event_start)
        print 'Processing summary information...'
        event_start = time.time()

        ## Create summary table
        print_summary(gen_shelf,symbol_table,'generation_%03d_summary.txt'%(n+1))

        event_end = time.time()
        gen_end = time.time()
        print
        print 'Time elapsed:',(event_end - event_start)
        print 'Time for this generation:',(gen_end - gen_start)
        print

        ## Promote to next generation
        last_gen = gen_shelf

    end_time = time.time()
    print 'Total elapsed time:',(end_time - init_time)

    for proc in procs:
        comm.send('EXIT',dest=proc,tag=2)

else:
    symbol_table = comm.bcast(symbol_table,root=0)

    if not symbol_table:
        sys.stdout.write('Process %d exiting.\n'%(mpi_rank))
        sys.exit()

    sys.stdout.write('Process %d runs.....\n'%(mpi_rank))
    
    ## Make a stack
    stack = deque()

    while 1:
        ## Clean up
        gc.collect()

        ## Send info on readiness
        comm.send(mpi_rank,dest=0,tag=1)

        ## Grab a message
        message = comm.recv(source=0,tag=2)

        if message == 'EXPAND':
            sys.stdout.write('Process %d expanding.....\n'%(mpi_rank))
            state = comm.recv(source=0,tag=3)
            base_prob_dict = comm.recv(source=0,tag=4)
        elif message == 'COMBINE':
            sys.stdout.write('Process %d combining.....\n'%(mpi_rank))
            comm.send(None,dest=0,tag=3)
        elif message == 'WAIT':
            sys.stdout.write('Process %d waiting.....\n'%(mpi_rank))
            message = comm.recv(source=0,tag=3)
        elif message == 'EXIT':
            break

        

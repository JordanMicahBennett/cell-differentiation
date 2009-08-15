# General discrete cell simulator program
# This file implements the Stem Cell model described in a classic paper
# "A Stochastic Model of Stem Cell Proliferation, Based on the Growth
# of Spleen Colony-Forming Cells"
# J.E. Till, E.A. McCulloch, and L. Siminovitch, PNAS 51:29-36 (1964)
# Program by Prof. M. Colvin, UC Merced
# This program is in the public domain for all purposes

# Bring in useful modules
from __future__ import division

import random
import math
import sys

# Set simulation parameters:
ntrials=1000         
nwells=1000          # Number of separate wells
ngen=8             # Number of generations
n_initial_cfcs=1    # Number CFC cells engrafted
p2=0.61             # CFC division probability

# Print parameter:
print "Number of Monte Carlo Trials:",nwells
print "Number of Generations:",ngen
print "Initial number of CFCs:",n_initial_cfcs
print "CFC division probability (P2):",p2

# The class cell holds info about a particular cell
class cell:
    def __init__(self, type):
        self.type=type
    def printcell(self): print self.type
    def __eq__(self, acell):
        return (self.type==acell.type)
            
# The class population holds the population of cells
class population:
    def __init__(self, types, fatedict):
        self.fatedict=fatedict
        self.types=types
        self.pop=[]
        self.ngen=0
    def printpopfates(self):
        print "Cells fates in population:"
        for type in self.types:
            print type
            self.fatedict[type].printfates()
    def printpop(self):
        print "Population:"
        i=1
        for cell in self.pop:
            print i, 
            cell.printcell()
            i+=1
    def addcell(self, cell): self.pop.append(cell)
    def update(self):
        oldpop=self.pop
        self.pop=[]
        for acell in oldpop:
            self.pop+=self.fatedict[acell.type].choosefate()

    def countcells(self,acell):
        icount=0
        for bcell in self.pop:
            if (acell.type==bcell.type): icount+=1 #Testing efficiency
#            if (acell == bcell): icount+=1
        return icount            
                
# The class fates describes all the ways a cell changes each generation
class fates:
    def __init__(self):
        self.fatelist=[]
    def addfate(self,fate):
        self.fatelist.append(fate)
    def nfates(self): print len(self.fatelist)
    def printfates(self):
        for i in self.fatelist: i.printfate()
    def choosefate(self):
        rand=random.random()
        i=0
        prob=self.fatelist[i].prob
        while(rand > prob):
            i+=1
            prob+=self.fatelist[i].prob
        #self.fatelist[i].printfate()
        newcells=[]
        for atype in self.fatelist[i].newtypes:
            newcells.append(cell(atype))
        return newcells
        
# The class fate describes one particular fate
class fate:
    def __init__(self,newtypes,prob):
        self.newtypes=newtypes
        self.prob=prob
    def printfate(self): print self.newtypes, self.prob

# The stats class keeps univariate stats on a variable
class stats:
    def __init__(self):
        self.n=0
        self.sumx=0.0
        self.sumx2=0.0
        self.sumx3=0.0
    def add(self, x):
        self.n+=1
        self.sumx+=x
        self.sumx2+=x*x
        self.sumx3+=x*x*x
    def ave(self):
        return float(self.sumx)/float(self.n)
    def stdev(self):
        return math.sqrt((self.sumx2-(self.sumx*self.sumx)/float(self.n))
                         /float(self.n-1))

# Utility class for building histograms
class histogram:         
    def __init__(self, nbins, binstart, binend):
        self.maxbars=30     #Maximum number bars in printed histogram
        self.nbins=nbins
        self.binstart=binstart
        self.binend=binend
        binwidth=(float(binend)-float(binstart))/float(nbins)
        self.bins=[]
        self.sum_x=0.0
        self.sum_x2=0.0
        self.n=0
        self.nlow=0
        self.nhigh=0
        self.bins.append(binstart)
        for i in range(nbins):
            self.bins.append((i+1)*binwidth)
        self.data=[0 for i in range(nbins)]
    def add(self,data):
        if (data < self.binstart):
            self.nlow+=1
            return
        if (data >= self.binend):
            self.nhigh+=1
            return
        i=1
        while (data >= self.bins[i]): i+=1
        self.data[i-1]+=1
        self.sum_x+=data
        self.sum_x2+=data*data
        self.n+=1
    def print_data(self, file=sys.stdout):
        file.write("n: %4d  nlow: %4d  nhigh: %4d\n" %
                   (self.n, self.nlow, self.nhigh))
        for i in range(self.nbins):
                if (i+1)%2: file.write("%2.2lf-%2.2lf: %-4d\n" %
                       (self.bins[i],self.bins[i+1], self.data[i]))
    def print_hist(self):
        maxdata=0;
        mean=float(self.sum_x)/float(self.n)
        stdev=math.sqrt(float(self.n)*self.sum_x2-self.sum_x*self.sum_x)/float(self.n)
        for i in range(self.nbins):
           if self.data[i]>maxdata:
               maxdata=self.data[i]
        print 'Mean: %4.2lf' % (mean)
        print 'Standard deviation: %4.2lf' % (stdev)
        for i in range(self.nbins):
            bars=int(self.maxbars*float(self.data[i])/float(maxdata))
            print '%2.2lf-%2.2lf: %-4d %s' % (self.bins[i],self.bins[i+1],
                                                   self.data[i],bars*"=")
# Utility class for building histograms
class histogram_cumulative:         
    def __init__(self, nbins, binstart, binend):
        self.maxbars=30     #Maximum number bars in printed histogram
        self.nbins=nbins
        self.binstart=binstart
        self.binend=binend
        binwidth=(float(binend)-float(binstart))/float(nbins)
        self.bins=[]
        self.sum_x=0.0
        self.sum_x2=0.0
        self.n=0
        self.nlow=0
        self.nhigh=0
        self.bins.append(binstart)
        for i in range(nbins):
            self.bins.append((i+1)*binwidth)
        self.data=[0 for i in range(nbins)]
    def add(self,data):
        if (data < self.binstart):
            self.nlow+=1
            return
        if (data >= self.binend):
            self.nhigh+=1
            return
        i=self.nbins
        while (data < self.bins[i]):
            self.data[i-1]+=1
            self.sum_x+=data
            self.sum_x2+=data*data
            self.n+=1
            i-=1
    def print_data(self, file=sys.stdout):
        file.write("n: %4d  nlow: %4d  nhigh: %4d\n" %
                   (self.n, self.nlow, self.nhigh))
        for i in range(self.nbins):
            if (i+1)%2: file.write("%2.2lf-%2.2lf: %-4d\n" %
                       (self.bins[i],self.bins[i+1], self.data[i]))
    def print_hist(self, file=sys.stdout):
        maxdata=0;
        mean=float(self.sum_x)/float(self.n)
        stdev=math.sqrt(float(self.n)*self.sum_x2-self.sum_x*self.sum_x)/float(self.n)
        for i in range(self.nbins):
           if self.data[i]>maxdata:
               maxdata=self.data[i]
        file.write("Mean: %4.2lf\n" % (mean))
        file.write("Standard deviation: %4.2lf\n" % (stdev))
        for i in range(self.nbins):
            bars=int(self.maxbars*float(self.data[i])/float(maxdata))
            file.write("%2.2lf-%2.2lf: %-4d %s\n" % (self.bins[i],self.bins[i+1],
                                           self.data[i],bars*"="))

n_zero_cfcs=[] # Number of runs that produce no cfcs
cfc_count=[]
for i in range(ngen):
    n_zero_cfcs.append(0)
    cfc_count.append(0)

# Defining the cells
# In this model we have 2 cell types:
# cfc: Colony Forming Cells
# edc: Early Differentiated cells

# This list of types will be used in keeping statistics
types=['cfc', 'edc']

# cfc fates:
fatecfc=fates()
# with probability p2, produce 2 cfc:
#p2=0.6
fatecfc.addfate(fate(['cfc','cfc'],p2))
# with probability p0, produce 1 edc:
p0=1-p2
fatecfc.addfate(fate(['edc'],p0))
#print "CFC Fates"
#fatecfc.printfates()

# edc fates:
fateedc=fates()
fateedc.addfate(fate(['edc'],1))
#print "EDC Fates"
#fateedc.printfates()

# Add this to the fate dictionary
fatedict={'cfc':fatecfc, 'edc':fateedc}

# keep some stats
cfc_hist=histogram(1000, 0, 1000)
#edc_hist=histogram(750, 0, 750)
#cfc_hist_cum=histogram_cumulative(750, 0, 750)

#Initialize stats classes
cfc_stats=stats()
edc_stats=stats()
cfc_gen_stats=[stats() for x in range(ngen)]
edc_gen_stats=[stats() for x in range(ngen)]

# Outer loop of number of simulations
for isim in range(nwells):
    ##print "Simulation trial:", isim
    if (isim%(nwells/10))==0:
        print "=",
    # First clear the stats and population
    pop=population(types, fatedict)
    #pop.printpopfates()
    for i in range(n_initial_cfcs):
        pop.addcell(cell('cfc'))

    # Now lets run the simulation for ngen generations
    for igen in range(ngen):
        # print "generation", igen
        # pop.printpop()
        # print "Number cfc cells;", pop.countcells(cell('cfc'))
        pop.update()
        ncfc=pop.countcells(cell('cfc'))
        nedc=pop.countcells(cell('edc'))
        cfc_gen_stats[igen].add(ncfc)
        edc_gen_stats[igen].add(nedc)
        # Increment counter if no cfc's present
        if (ncfc==0):
            n_zero_cfcs[igen]+=1
        cfc_count[igen]=ncfc

#    for igen in range(ngen):
#        print '%d'%(cfc_count[igen]),
#    print
    
# Update statistics and histogram classes
    cfc_stats.add(ncfc)
    edc_stats.add(nedc)
    cfc_hist.add(ncfc)
#    edc_hist.add(nedc)
#    cfc_hist_cum.add(ncfc)
    
# Print out relevant stats
print "\nCFC Average:",cfc_stats.ave(),"SD:", cfc_stats.stdev()
print "EDC Average:",edc_stats.ave(),"SD:", edc_stats.stdev()
print "\nPercent not surviving (no cfc's) at each generation:"
for i in range(ngen):
    print "Gen=",i+1,"% not surviving:",(n_zero_cfcs[i]/nwells)

# Make generational tables    
print "\nGenerational ave results"
print "Gen","Ave CFC","SD CFC","Ave EDC","SD EDC"
for i in range(ngen):
    print '%d %6.2f %6.2f %6.2f %6.2f' % (i+1,cfc_gen_stats[i].ave(),cfc_gen_stats[i].stdev(),
                                          edc_gen_stats[i].ave(),edc_gen_stats[i].stdev())

# Print Histograms
#outfile=open("histograms.txt","w")
#cfc_hist.print_data(sys.stdout)
#edc_hist.print_data(outfile)
#outfile.write("CFC Cumulative Histogram\n")
#cfc_hist_cum.print_data(outfile)
#outfile.close()

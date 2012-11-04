#!/usr/bin/env python

# qmc: a naive implementation of the Quine-McCluskey algorithm.
# The algorithm optimizes two level multi-output boolean functions.
# Input: function's ONset (the set of minterms)
#                   and DCset (don't care set)
# Output: boolean form of the optimized function(s)
# Bibliography for the algorithm: C.Bolchini et al, 'Reti Logiche', Apogeo
#
# Copyright (C) 2011 Marcello Pogliani
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import math

VERBOSE = False # enables various prints during the algorithm execution. Used for debugging purpose or for the purpose of simulating the algorithm by hand and checking the various steps against those made by this program.
LITERALS_COST_FUNCTION = 0 # TODO explain (refer to Reti Logice book for explaination) ;)
IMPLICANTS_COST_FUNCTION = 1
DEFAULT_COST_FUNCTION = LITERALS_COST_FUNCTION

class Implicant:
    def __init__(self, element, num_of_vars, cost_function, coverset=False, dcmask=0, funmask=1):
        self.onmask = element
        self.dcmask = dcmask
        self.num_of_vars = num_of_vars
        self.funmask = funmask
        self.reduced = False
        self.cost_function = cost_function
        self.already_selected = False
        self.name = None
        if not coverset:
            self.coverset = {element}
        else:
            self.coverset = coverset

    def onesNumber(self):
        """ Number of 'ones' in the implicant (seen as a binary string) """
        onemask = self.onmask & ~self.dcmask
        ones_number = 0
        while onemask > 0:
            ones_number += onemask & 1
            onemask >>= 1
        return ones_number
    
    def defaultCost(self):
        """ The basic (without considering that the implicant could be shared)
        cost of the implicant, computed with the cost function set above """
        if self.cost_function == LITERALS_COST_FUNCTION:
            cost = self.num_of_vars
            for i in range(cost):
                if (self.dcmask >> i) & 1 == 1:
                    cost -= 1
            return cost
        else:
            return 1    # each implicant costs 1

    def setUpdatedCost(self):
        """ Update the cost, bringing it to that of shared implicants """
        self.already_selected = True

    def nextCost(self):
        """ cost of selecting another time the implicant """
        if self.already_selected:
            return self.sharedCost()
        else:
            return self.defaultCost()

    def sharedCost(self):
        if self.cost_function == LITERALS_COST_FUNCTION:    
            return 1
        else:
            return 0

    def joinedWith(self, implicant2):
        """ Create a new implicant by joining two implicants
        WARNING: this method doesn't check if they're joinable!!! """
        return Implicant(element = self.onmask & implicant2.onmask,
                         num_of_vars = max(self.num_of_vars, implicant2.num_of_vars),
                         cost_function = self.cost_function,
                         dcmask = self.dcmask | (self.onmask ^ implicant2.onmask),
                         coverset = self.coverset | implicant2.coverset,
                         funmask = self.funmask & implicant2.funmask )

    def isInListOfImplicants(self, list_of_implicants):
        """ Check whether the implicant is in a list. Very bad & ugly :-) """
        for element in list_of_implicants:
            if element.onmask == self.onmask and element.dcmask == self.dcmask:
                return True
        return False

    def hammingDistanceFrom(self, imp2):
        """ Computes the Hamming distance between two implicants,
        taking into account the don't cares """
        diffmask = self.onmask ^ imp2.onmask & ~(self.dcmask | imp2.dcmask)
        # values in either the 1st or the 2nd dcmask taken into account later
        distance = 0
        while diffmask > 0:
            distance += diffmask & 1
            diffmask >>= 1
        if self.dcmask != imp2.dcmask:  # difference in the dcmask is relevant
            diffmask = self.dcmask | imp2.dcmask
            while diffmask > 0:
                distance += diffmask & 1
                diffmask >>= 1
        return distance
        
    def maskToString(self, num_of_bits = None):
        """ A string, with the function mask in binary form (DRAFT...)
        WARNING: * zeroes on the right (after the last one) are MISSING! """
        t_str, k = '', 0
        num = self.funmask
        while num > 0:
            if num & 1 == 1:
                t_str = t_str + '1'
            else:
                t_str = t_str + '0'
            num >>= 1
            k += 1
        # pad with zeroes
        if num_of_bits != None and num_of_bits > k:
            t_str += '0' * (num_of_bits - k)
        return t_str

    def __str__(self):
        """ A string with the implicant in boolean form """
        t_str = ''
        for i in range(self.num_of_vars):
            if (self.dcmask >> self.num_of_vars - i - 1) & 1 == 0:
                ltr = chr(ord('a') + i)
                if (self.onmask >> self.num_of_vars - i - 1) & 1 == 1:
                    t_str += ltr
                else:
                    t_str += ltr + '\''
        return t_str
    
    # the following are just getters and\or setters, implemented 
    # only for the sake of making the rest of the code clearer
    def setReduced(self, isTrue = True):
        self.reduced = isTrue

    def isReduced(self):
        return self.reduced

class QmcFunction:
    def __init__(self, ONset, DCset):
        self.var_num = self.variableSetCardinality(ONset, DCset)
        self.ONset = ONset
        self.DCset = DCset
        self.to_be_covered = ONset

    def assignNumber(self, number):
        # will throw a LOT of exceptions if number is not assigned!!!!!!!!
        self.number = number

    def reset(self):
        """ revert the minterms to cover to the complete ONset """
        self.to_be_covered = self.ONset

    def variableSetCardinality(self, ONset, DCset):
        """ number of variables of the function """
        return int(math.ceil(math.log(max(ONset + DCset)+1, 2)))

    def updateMintermsToCover(self, implicant):
        """ Select a new implicant for coverage """
        self.to_be_covered = filter(lambda x : x not in implicant.coverset, self.to_be_covered)

class QuineMcCluskey:
    def __init__(self, function_list, cost_function):
        # some assertions...
        assert cost_function == IMPLICANTS_COST_FUNCTION \
           or cost_function == LITERALS_COST_FUNCTION
        # initializations & structures...
        self.cost_function = cost_function
        for k, fun in enumerate(function_list):
            fun.assignNumber(k)
        self.functions_to_simplify = function_list
        self.implicant_list = self.buildMintermList()
        self.sol = Solution(len(self.functions_to_simplify))
        # sort the list wrt the number of ones
        # (the list is assumed sorted by some loops later!!!)
        self.implicant_list.sort(key = lambda e : e.onesNumber())

    def numberOfVariables(self):
        n = 0
        for fun in self.functions_to_simplify:
            if fun.var_num > n:
                n = fun.var_num
        return n

    def buildFunMask(self, functions_covering):
        mask = 0
        for element in functions_covering:
            tmp = element
            mask |= 1 << tmp
        return mask

    def buildMintermList(self):
        """ Build the list of minterms that will be expanded given the ONset
        and DCset as lists of number representing in decimal notation the
        input values for which the output should be 1 (or don't cares) """
        nov = self.numberOfVariables()
        masks = {} # {element of ON\DCset:[list of functions it covers]}
        for fun in self.functions_to_simplify:
            for el in fun.ONset + fun.DCset:
                try:
                    masks[el].append(fun.number)
                except KeyError:
                    masks[el] = [fun.number]
        return [Implicant(key, nov, 
                          cost_function = self.cost_function,
                          funmask=self.buildFunMask(masks[key])) 
                for key in masks]

    def expansionStep(self, implicant_list):
        """ single step of the expansion procedure """
        newimp = []
        for numi1, i1 in enumerate(implicant_list):
            onesnumber = i1.onesNumber()
            for i2 in implicant_list[numi1 + 1:]:
                if i2.onesNumber() > onesnumber + 1:  # if I have two more ones
                # the hamming distance is >= 1 if the implicant list IS ORDERED
                    break
                joinfunmask = i1.funmask & i2.funmask
                if i1.hammingDistanceFrom(i2) == 1 and joinfunmask != 0:
                    if i1.funmask == joinfunmask:
                        i1.setReduced()
                    if i2.funmask == joinfunmask:
                        i2.setReduced()
                    tmp = i1.joinedWith(i2)
                    if not tmp.isInListOfImplicants(newimp):
                        newimp.append(tmp)
        return newimp

    def findPrimeImplicants(self):
        """ performs the expansion procedure and finds the prime implicants
        Writes the result in self.implicant_list (not to waste bytes...) """
        newimp = self.implicant_list
        while newimp != []:
            newimp = self.expansionStep(newimp)
            self.implicant_list += newimp
        self.implicant_list[:] = [i for i in self.implicant_list if not i.isReduced()]
        self.purgeTable()   # remove implicants covering only don't cares
        self.giveNameToImplicants()
        return self.implicant_list  # those are the prime implicants

    def giveNameToImplicants(self):
        n = 0
        for imp in self.implicant_list:
            imp.name = chr(ord('A') + n)
            n += 1

    def essentialityStep(self):
        """ A single essentiality step. Returns the set of covered
        minterms and the set of essential implicants that have been found. """
        done_something = False
        for fun in self.functions_to_simplify:
            selected = self.doEssentialityForFunction(fun)
            for implicant in selected:
                done_something = True
                # add the implicant to the solution & related things
                implicant.setUpdatedCost()
                fun.updateMintermsToCover(implicant)
                self.sol.addImplicant(implicant, fun)
        # if the implicants we've selected can't cover anything more
        # they should be removed from the table
        self.purgeTable()
        return done_something

    def doEssentialityForFunction(self, fun):
        """ The essentiality step for a single function """
        selected = set()
        for column in fun.to_be_covered:
            selected_implicant = self.implicantEssentialForElementInColumn(column, fun)
            if selected_implicant:
                selected.add(selected_implicant)
        return selected

    def implicantEssentialForElementInColumn(self, column, fun):
        """ Is there a prime implicant essential for the minterm 'column' 
        in the function 'fun_num'? If so, returns it """
        for number, implicant in enumerate(self.implicant_list):
            if column in implicant.coverset \
               and (implicant.funmask >> fun.number) & 1 == 1:
                for implicant2 in self.implicant_list[number + 1:]:
                    if column in implicant2.coverset \
                       and (implicant2.funmask >> fun.number) & 1 == 1:
                        return False
                return implicant

    def purgeTable(self):
        """ Purge the table (implicant_list) from empty rows """
        todelete = set()
        for implicant in self.implicant_list:
            todelete.add(implicant)
            for f in self.functions_to_simplify:
                if (implicant.funmask >> f.number) & 1 == 1:
                    s = {i for i in f.to_be_covered if i in implicant.coverset}
                    if s != set():
                        todelete.remove(implicant)
                        break
        self.implicant_list = [i for i in self.implicant_list if i not in todelete]

    def rowDominanceStep(self):
        """ Purge the implicant list from dominated rows """
        todelete = set()
        done_something = False
        for i, implicant in enumerate(self.implicant_list):
            for implicant2 in self.implicant_list[i + 1:]:
                a_cbd = True
                b_cbd = True
                for f in self.functions_to_simplify:
                    a = set()
                    b = set()
                    if (implicant.funmask >> f.number) & 1 == 1:
                        a = {x for x in implicant.coverset 
                               if x in f.to_be_covered}
                    if (implicant2.funmask >> f.number) & 1 == 1:
                        b = {x for x in implicant2.coverset 
                               if x in f.to_be_covered}
                    if not b.issubset(a): a_cbd = False
                    if not a.issubset(b): b_cbd = False
                if a_cbd and implicant.nextCost() <= implicant2.nextCost():
                    todelete.add(implicant2)
                    done_something = True
                elif b_cbd and implicant2.nextCost() <= implicant.nextCost():
                    todelete.add(implicant)
                    done_something = True
        self.implicant_list = [i for i in self.implicant_list 
                                 if not i in todelete]
        return done_something

    def columnDominanceStep(self):
        """ The list of minterms to cover purged from dominated columns """
        done_something = False
        for f in self.functions_to_simplify:
            to_be_deleted = set()
            for numcol, col1 in enumerate(f.to_be_covered):
                for col2 in f.to_be_covered[numcol + 1:]:
                    lc1 = {imp for imp in self.implicant_list
                               if col1 in imp.coverset}
                    lc2 = {imp for imp in self.implicant_list
                               if col2 in imp.coverset}
                    if lc1.issubset(lc2): # col1 dominates col2
                        to_be_deleted.add(col2)
                        done_something = True
                    elif lc2.issubset(lc1): # col2 dominates col1
                        to_be_deleted.add(col1)
                        done_something = True
            f.to_be_covered = [i for i in f.to_be_covered 
                                 if i not in to_be_deleted]
        return done_something
    
    def branchAndBound(self):
        """ Branch And Bound algorithm to find the optimal coverage """
        set_s = [Solution(len(self.functions_to_simplify))]
        # set (list for simplicity) of partially explored solutions
        solution = Solution(len(self.functions_to_simplify)) # empty solution
        bound = 'inf'       # means infinity
        while set_s:        # while we have solutions to explore
            k = set_s[0]    # select a partial solution from set_s (the first)
            set_s = set_s[1:] # delete set_s[0] from the solution set
            if bound == 'inf' or k.getCost() < bound:
                branched_list = self.branch(k)
                for k_i in branched_list:
                    b = k_i.getCost()
                    if bound == 'inf' or b < bound:
                        if k_i.isComplete(self.functions_to_simplify):
                            solution = k_i
                            bound = b
                            if VERBOSE:
                                print ' |=> Complete solution. New bound = ', b
                        else:
                            set_s.append(k_i)
        return solution

    def branch(self, solution):
        """ Branches a _partial_ solution by selecting another implicant """
        k_list = []
        for f in self.functions_to_simplify:
            for i in f.to_be_covered:
                if not solution.covers(i, f):
                    for j in self.implicant_list:
                        if i in j.coverset and (j.funmask >> f.number) & 1 == 1:
                            k_list.append(solution.branchWith(j, f.number))
                    return k_list

    def simplify(self):
        """ Start all the optimization process! """
        if not self.reduceTable():
            if VERBOSE:
                print 'Branch and bound'
            solution = self.branchAndBound()
            self.sol.mergeWith(solution)
        return self.sol

    def reduceTable(self):
        """ Try to simplify the problem by applying essentiality & dominance
        criteria before using the branch and bound algorithm """
        done_something = True
        while done_something:
            if VERBOSE: print 'Detecting essential implicants:'
            done_something = self.essentialityStep()
            
            if done_something:
                found_complete_cover = self.sol.isComplete(self.functions_to_simplify)
                
                if VERBOSE and not found_complete_cover:
                    self.printTable()

            elif VERBOSE: print '(essential implicants not found)'
            if VERBOSE: print 'Simplifying table by deleting dominated rows & columns'
            
            tmp1 = self.rowDominanceStep()
            tmp2 = self.columnDominanceStep()
            done_something = tmp1 or tmp2 or done_something

            if VERBOSE:
                if tmp1 or tmp2:
                    print 'Resulting table:'
                    self.printTable()
                if not tmp1: print '(dominated rows not found)'
                elif not tmp2: print '(dominated columns not found)'
        return found_complete_cover

    def printImplicantList(self):
        """ prints the specified list of implicants """
        for implicant in self.implicant_list:
            if implicant.isReduced():
                selected_string = 'X'
            else:
                selected_string = ' '
            print '%3c = %s %c\t' % (implicant.name, 
                            implicant, selected_string),
            print implicant.maskToString(len(self.functions_to_simplify)), '\t',
            for i in implicant.coverset:
                print i,
            print '' # newline

    def printTable(self):
        """ prints the coverage table in a decent format. C column = cost"""
        print ' ' * 8,
        for f in self.functions_to_simplify:
            for column in f.to_be_covered:
                print "%3d  " % column,
            print '|',
        print ' C'
        for implicant in self.implicant_list:
            print '%3c ->' % implicant.name,
            for f in self.functions_to_simplify:
                print '|',
                for column in f.to_be_covered:
                    if column in implicant.coverset \
                       and (implicant.funmask >> f.number) & 1 == 1:
                            print '  X  ',
                    else:
                        print '     ',
            print '| ',implicant.nextCost()

class Solution():
    def __init__(self, fnum = 0):
        self.selected = []
        self.cost = 0
        for k in range(fnum):
            self.selected.append(set())

    def getCost(self):
        return self.cost
    
    def recomputeCost(self):
        """ computes the cost of the partial solution """
        cost = 0
        sel = set()
        for set_f in self.selected:
            for implicant in set_f:
                if implicant in sel:    # shared implicant
                    if implicant.cost_function == LITERALS_COST_FUNCTION:
                        cost += 1
                else:
                    sel.add(implicant)
                    cost += implicant.defaultCost()
        return cost

    def covers(self, target_imp, f):
        """ Is this solution covering target_imp wrt the fnum function? """
        for implicant in self.selected[f.number]:
            if target_imp in implicant.coverset:
                return True
        return False

    def isComplete(self, functions_to_simplify):
        """ Is this solution a complete one? """
        for f in functions_to_simplify:
            fcv = set()
            for im in self.selected[f.number]:
                fcv |= im.coverset
            for n in f.to_be_covered:
                if n not in fcv:
                    return False
        return True

    def branchWith(self, implicant, fnum):
        """ Create a new Solution object by adding implicant for fnum """
        t = Solution()
        for x in self.selected:
            t.selected.append(set() | x)
        t.cost = self.cost
        t.selected[fnum].add(implicant)
        # update the cost...
        updated_cost = False
        if implicant.nextCost() == implicant.defaultCost():
            # check if the implicant is already selected for some other functions...
            for fun in range(len(self.selected)):
                if fun != fnum and implicant in self.selected[fun]:
                    t.cost += implicant.sharedCost()
                    updated_cost = True
                    break
        if not updated_cost:
            t.cost += implicant.nextCost()
        return t
    
    def mergeWith(self, sol2):
        for i in range(len(self.selected)):
            self.selected[i] |= sol2.selected[i]
        # this code is executed only one time
        # so it's quicker just to recompute the cost from the beginning
        self.cost = self.recomputeCost()
    
    def addImplicant(self, implicant, fun):
        """ add an implicant to this Solution object """
        self.selected[fun.number].add(implicant)
        self.cost = self.recomputeCost()
        
    def __str__(self):
        """ returns a string containing the SOP form of the logical function
        of all the implicants selected until now or specified in cov """
        a, b, t_str = [], [], ''
        for num, selected_for_cur_f in enumerate(self.selected):
            t_str += 'F%d = ' % (num + 1)
            is_first_step = True
            for element in selected_for_cur_f:
                if not is_first_step:
                    t_str += ' + '
                else: is_first_step = False
                
                if element in a:
                    t_str += element.name
                    b.append(element)
                else:
                    a.append(element)
                    found = False
                    for selected_for_next_f in self.selected[num+1:]:
                        if element in selected_for_next_f:
                            found = True
                            t_str += element.name
                    if not found:
                        t_str += element.__str__()
            t_str += '\n'
        b.sort(key = lambda x: x.name)
        for element in b:
            t_str += element.name + ' = ' + element.__str__() + '\n'
        return t_str

# TODO parse command line args as input instead than modifying the source code!
if __name__ == '__main__':
    # an example of use of the library
    VERBOSE = True # enable debug prints when running with the example...
    q = QuineMcCluskey([
                       QmcFunction([3,4,6,7], [0,2]),
                       QmcFunction([3,4,6], [2,5]),
                       QmcFunction([0,2], [4,6])
                       ], LITERALS_COST_FUNCTION)
    print ' --- [Expansion] Prime implicants: ---'
    q.findPrimeImplicants()
    q.printImplicantList()
    q.printTable()
    print 'Solution found:'
    print q.simplify()
    print 'Cost:', q.sol.getCost()

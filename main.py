from bessel import truncatedBessel, dynamicBessel, prettyErrorTable
from bisection import bisectionMethod
from secant import secantMethod
from newton import newtonsMethod
from tableprinter.tableprinter import prettyTable
from scipy.special import jv
import sys

#CONSTANTS
## Our good friend delta:
delta = 10**(-16)
## For calculating J_0:
xVals = [1,5,10,20]
nVals = [10,25,50,75,100]
## For finding roots of J_0:
intervals = [[1,3],[4,6],[8,10]]
newtonGuesses = [3,6,10]
epsilons = [10**(-x) for x in range(6,11)]

def getResultTable(method, intervals, epsilons):
    besselZero = lambda x: jv(0,x)
    specifiedMethod = lambda x,y: method(besselZero, x, 200, y, delta)
    # Yes, this next bit of code looks pretty bad, but it generates perfectly formatted tables.
    # I tried refactoring it, but honestly it's pretty straightforward already, just long.
    return [   ["[{0},{1}]".format(domain[0],domain[1])] +
               ["{0:.13}, {1} iters.".format(specifiedMethod(domain,ep)[0], specifiedMethod(domain,ep)[2]) for ep in epsilons]
           for domain in intervals]

def generatePrettyTable(method, intervals, epsilons):
    return prettyTable(getResultTable(method, intervals, epsilons), epsilons, "[a,b]", "eps.")

# Newton's method is special and needs its own method for making a result table.
def getNewtonTable(guesses, epsilons):
    besselZero = lambda x: jv(0,x)
    besselPrime  = lambda x: -jv(1,x)
    boundNewton = lambda x,y: newtonsMethod(besselZero, besselPrime, x, 200, y, delta)
    return [[str(guess)] + ["{0:.13}, {1} iters.".format(boundNewton(guess, ep)[0],boundNewton(guess,ep)[2]) for ep in epsilons] for guess in guesses]

def generateNewtonTable(guesses, epsilons):
    return prettyTable(getNewtonTable(guesses, epsilons), epsilons, "Guess", "eps.")

def returnTables():
    print("Table of relative errors for the truncated Bessel function:")
    print(prettyErrorTable(xVals, nVals, truncatedBessel))
    print("\n")
    print("Table of relative errors for the recursive Bessel function:")
    print(prettyErrorTable(xVals, nVals, dynamicBessel))
    print("\n"*3)
    print("Bisection method results for J_0:")
    print(generatePrettyTable(bisectionMethod, intervals, epsilons))
    print("\n")
    print("Secant method results for J_0:")
    print(generatePrettyTable(secantMethod, intervals, epsilons))
    print("\n")
    print("Newton's method results for J_0:")
    print(generateNewtonTable(newtonGuesses, epsilons))

def main():
    """ I'm having this just print the tables for now, but I'm leaving myself the option to
        add functionality in the future. """
    returnTables()

if __name__ == '__main__':
    main()

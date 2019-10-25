from bessel import truncatedBessel, dynamicBessel
from bisection import bisectionMethod
from secant import secantMethod
from newton import newtonsMethod
from tableprinter.tableprinter import prettyTable
from scipy.special import jv
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

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

## BESSEL FUNCTION ERROR CALCULATIONS

def  relErr(x, n, callback):
    """ Returns the relative error for the calculation of the
        Bessel function of order zero at the point x using the function passed
        as the callback. """
    correctValue = jv(0,x)
    estimate = callback(x,n)
    return abs(correctValue - estimate)/abs(correctValue)

def calculateErrorTable(xvals, nvals, method):
    """ This looks like a complicated one-line expression, but really all it does
        is construct a 2-dimensional matrix which has rows where the first element
        is the value of x we are testing and the remaining elements are the relative
        errors for those tests at various values of N, i.e. just a table of values."""
    return [[x] + [relErr(x,n,method) for n in nvals] for x in xvals]

def prettyErrorTable(xVals, nVals, method):
    tableMatrix = calculateErrorTable(xVals, nVals, method)
    return prettyTable(tableMatrix, nVals, "x", "N")

## ROOT FINDING METHODS

def getResultTable(method, intervals, epsilons):
    besselZero = lambda x: jv(0,x)
    specifiedMethod = lambda x,y: method(besselZero, x, 200, y, delta)
    # Yes, this next bit of code looks pretty bad, but it generates perfectly formatted tables.
    # I tried refactoring it, but honestly it's pretty straightforward already, just long.
    return [   ["[{0},{1}]".format(domain[0],domain[1])] +
               ["{0:.13}, {1} iters.".format(specifiedMethod(domain,ep)[0], specifiedMethod(domain,ep)[2]) for ep in epsilons]
           for domain in intervals]

def getResultTableForPlot(method, intervals, epsilons):
    besselZero = lambda x: jv(0,x)
    specifiedMethod = lambda x,y: method(besselZero, x, 200, y, delta)
    return [[specifiedMethod(domain, ep)[2] for ep in epsilons] for domain in intervals]

def generatePrettyTable(method, intervals, epsilons):
    return prettyTable(getResultTable(method, intervals, epsilons), epsilons, "[a,b]", "eps.")

# Newton's method is special and needs its own method for making a result table.
def getNewtonTable(guesses, epsilons):
    besselZero = lambda x: jv(0,x)
    besselPrime  = lambda x: -jv(1,x)
    boundNewton = lambda x,y: newtonsMethod(besselZero, besselPrime, x, 200, y, delta)
    return [[str(guess)] + ["{0:.13}, {1} iters.".format(boundNewton(guess, ep)[0],boundNewton(guess,ep)[2]) for ep in epsilons] for guess in guesses]

def getNewtonTableForPlot(guesses, epsilons):
    besselZero = lambda x: jv(0,x)
    besselPrime  = lambda x: -jv(1,x)
    boundNewton = lambda x,y: newtonsMethod(besselZero, besselPrime, x, 200, y, delta)
    return [[boundNewton(guess, ep)[2] for ep in epsilons] for guess in guesses]

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

def exportFigures(filepath):
    filename, extension = os.path.splitext(filepath)
    truncatedError = calculateErrorTable(xVals, nVals, truncatedBessel)
    dynamicError = calculateErrorTable(xVals, nVals, dynamicBessel)
    combined = [[truncatedError[n][1:]] + [dynamicError[n][1:]] for n in range(4)]
    plt.figure()

    for k in range(4):
        plt.subplot(221 + k)
        plt.plot(nVals, combined[k][0], "b")
        plt.plot(nVals, combined[k][1], "r")
        plt.title("$x = {}$".format(xVals[k]))
        plt.xlabel("Number of terms")
        plt.ylabel("Relative error")
        plt.yscale('log')

    plt.subplots_adjust(wspace = 0.5, hspace=0.7)

    plt.savefig(filename + "_bessel" + extension)

    plt.close()

    bisectionPlotTable = getResultTableForPlot(bisectionMethod, intervals, epsilons)
    secantPlotTable = getResultTableForPlot(secantMethod, intervals, epsilons)
    tables = [bisectionPlotTable, secantPlotTable]
    plt.figure()

    for k in range(len(intervals)):
        plt.subplot(221 + k)
        plt.plot(epsilons, tables[0][k], "b")
        plt.plot(epsilons, tables[1][k], "r")
        plt.title("Interval: " + str(intervals[k]))
        plt.xlabel("Epsilon")
        plt.ylabel("Iterations")
        plt.xscale('log')

    plt.subplots_adjust(wspace=0.6, hspace=0.7)
    plt.savefig(filename + "_bisection_secant" + extension)
    plt.close()

    newtonTable = getNewtonTableForPlot(newtonGuesses, epsilons)
    plt.figure()
    for k in range(len(newtonGuesses)):
        plt.subplot(221 + k)
        plt.plot(epsilons, newtonTable[k], "g")
        plt.title("Initial guess: " + str(newtonGuesses[k]))
        plt.xlabel("Epsilon")
        plt.ylabel("Iterations")
        plt.xscale("log")

    plt.subplots_adjust(wspace=0.6, hspace=0.7)
    plt.savefig(filename + "_newton" + extension)
    plt.close()

def usage():
    print("""
    Usage: python main.py [filepath.png]\n
    Add a filepath to export plots to that location.""")

def main():
    """ Arguments choose whether to print tables or export plots. """
    if len(sys.argv) == 1:
        returnTables()
    elif len(sys.argv) == 2:
        exportFigures(sys.argv[1])
    else:
        usage()


if __name__ == '__main__':
    main()

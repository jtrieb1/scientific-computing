from interpolators import makeEquiDistRange, makeCosRange, polynomialInterpolate
from spline import spline
import numpy as np
import matplotlib.pyplot as plt
import math
from tableprinter.tableprinter import prettyTable
from integration import compositeTrapezoidal, compositeSimpson, gaussianQuad
import sys

# Defining our constants
cVals = [1,4,25]
nVals = [6,10,14,18]

def function(c):
  assert c > 0
  return lambda x: 1 / (c*(x**2) + 1)

funcList = [function(c) for c in cVals]

def error(func1, func2, x):
  try:
    return abs(func1(x)-func2(x))
  except:
    return 0

def composedError(func1, func2):
  return lambda x: error(func1, func2, x)

def chartGen1(save, path=""):
  x1 = np.arange(-1.0,1.0,0.01)
  rows = len(funcList)
  colors = "bgrcmy"
  fig=plt.figure()
  base = rows * 100 + 21
  for i in range(rows):
    plt.subplot(base + 2*i)
    plt.title("Equidistant Point Interpolation")
    approxList = []
    plt.plot(x1, funcList[i](x1), "k")
    for n in nVals:
      approxList.append(polynomialInterpolate(funcList[i],makeEquiDistRange(n)))
    for j in range(len(approxList)):
      plt.plot(x1, approxList[j](x1), colors[j], label=str(nVals[j]))
    plt.subplot(base + 2*i + 1)
    plt.title("Error")
    for k in range(len(approxList)):
      plt.plot(x1, composedError(funcList[i],approxList[k])(x1), colors[k], label=str(nVals[k]))
  plt.subplots_adjust(wspace=1.,hspace=0.7)
  lines_labels = [fig.axes[0].get_legend_handles_labels()]
  lines, labels = [sum(lol,[]) for lol in zip(*lines_labels)]
  fig.legend(lines, labels, title="Intervals", loc="center")
  if save:
    plt.savefig(path + "-Equi.jpg")
  else:
    plt.show()

  fig=plt.figure()
  for i in range(rows):
    plt.subplot(base + 2*i)
    plt.title("Cosine-based Node Interpolation")
    approxList = []
    plt.plot(x1, funcList[i](x1), "k")
    for n in nVals:
      approxList.append(polynomialInterpolate(funcList[i],makeCosRange(n)))
    for j in range(len(approxList)):
      plt.plot(x1, approxList[j](x1), colors[j], label=str(nVals[j]))
    plt.subplot(base + 2*i + 1)
    plt.title("Error")
    for k in range(len(approxList)):
      plt.plot(x1, composedError(funcList[i],approxList[k])(x1), colors[k], label=str(nVals[k]))
  plt.subplots_adjust(wspace=1.,hspace=0.7)
  lines_labels = [fig.axes[0].get_legend_handles_labels()]
  lines, labels = [sum(lol,[]) for lol in zip(*lines_labels)]
  fig.legend(lines, labels, title="Intervals", loc="center")
  if save:
    plt.savefig(path + "-Cos.jpg")
  else:
    plt.show()

  fig = plt.figure()
  
  base = rows * 100 + 21
  for i in range(rows):
    plt.subplot(base + 2*i)
    plt.title("Spline Interpolation")
    plt.plot(x1,funcList[i](x1),"k")
    spList = [spline(funcList[i], makeEquiDistRange(n)) for n in nVals]
    for j in range(len(spList)):
      plt.plot(x1, list(map(spList[j],(x1))), colors[j], label=str(nVals[j]))
    plt.subplot(base + 2*i + 1)
    plt.title("Error")
    for k in range(len(nVals)):
      iSpline = spline(funcList[i], makeEquiDistRange(nVals[k]))
      spError = composedError(funcList[i], iSpline)
      plt.plot(x1, list(map(spError,(x1))), colors[k], label=str(nVals[k]))
  plt.subplots_adjust(wspace=1.,hspace=0.7)
  lines_labels = [fig.axes[0].get_legend_handles_labels()]
  lines, labels = [sum(lol,[]) for lol in zip(*lines_labels)]
  fig.legend(lines, labels, title="Intervals", loc="center")
  if save:
    plt.savefig(path + "-Spline.jpg")
  else:
    plt.show()

### Part 2 of project begins here.
nVals2 = [5,10,15,20,25,30]

def intError(callback, method, n, interval):
  return abs(callback(method, n, interval)-math.pi)

def errorTable(nVals, cbList, method, interval):
  return [[n] + ["{0:.10}, err: {1:.10}".format(callback(method, n, interval), intError(callback,method,n,interval)) for callback in cbList] for n in nVals]

def makePrettyTable(cbList, method, nVals, interval, names=["Trapezoidal","Simpson"]):
  return prettyTable(errorTable(nVals, cbList, method, interval), names, "N","method",padding=2)

cbList = [compositeTrapezoidal, compositeSimpson]

gnVals = [2,4,8,12,16,32]

def integral2(x):
  return 4 / (1 + x**2)

def integral1(x):
  return math.sqrt(4-(x**2))

def chartGen2():
  print(makePrettyTable(cbList, integral1, nVals2, (0,2)))
  print()
  print(makePrettyTable(cbList,integral2,nVals2,(0,1)))
  print()
  print()
  print(makePrettyTable([gaussianQuad],integral1,gnVals,(0,2),names=["Gaussian"]))
  print()
  print(makePrettyTable([gaussianQuad], integral2,gnVals,(0,1),names=["Gaussian"]))

def usage():
  print("""
  USAGE: python main.py [1|2] [PathToFolder].
  Given a number, returns the solutions to that section of the project. 
  If a path is provided, pictures of the generated graphs will be saved.""")

def main():
  argv = sys.argv
  argc = len(argv)
  if argc >= 4 or argc < 2:
    return usage()
  if argv[1] == "1" and argc == 3:
    return chartGen1(True, argv[2])
  elif argv[1] == "1":
    return chartGen1(False)
  elif argv[1] == "2":
    return chartGen2()

if __name__ == "__main__":
  main()
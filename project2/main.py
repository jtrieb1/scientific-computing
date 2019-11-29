from interpolators import makeEquiDistRange, makeCosRange, polynomialInterpolate
from spline import spline
import numpy as np
import matplotlib.pyplot as plt
import math
from tableprinter.tableprinter import prettyTable
from integration import compositeTrapezoidal, compositeSimpson, gaussianQuad

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

x1 = np.arange(-1.0,1.0,0.01)
rows = len(funcList)
colors = "bgrcmy"
plt.figure()
base = rows * 100 + 21
for i in range(rows):
  plt.subplot(base + 2*i)
  approxList = []
  plt.plot(x1, funcList[i](x1), "k")
  for n in nVals:
    approxList.append(polynomialInterpolate(funcList[i],makeEquiDistRange(n)))
  for j in range(len(approxList)):
    plt.plot(x1, approxList[j](x1), colors[j])
  plt.subplot(base + 2*i + 1)
  for k in range(len(approxList)):
    plt.plot(x1, composedError(funcList[i],approxList[k])(x1), colors[k])

plt.show()

plt.figure()
for i in range(rows):
  plt.subplot(base + 2*i)
  approxList = []
  plt.plot(x1, funcList[i](x1), "k")
  for n in nVals:
    approxList.append(polynomialInterpolate(funcList[i],makeCosRange(n)))
  for j in range(len(approxList)):
    plt.plot(x1, approxList[j](x1), colors[j])
  plt.subplot(base + 2*i + 1)
  for k in range(len(approxList)):
    plt.plot(x1, composedError(funcList[i],approxList[k])(x1), colors[k])

plt.show()

plt.figure()
base = rows * 100 + 21
for i in range(rows):
  plt.subplot(base + 2*i)
  plt.plot(x1,funcList[i](x1),"k")
  spList = [spline(funcList[i], makeEquiDistRange(n)) for n in nVals]
  for j in range(len(spList)):
    plt.plot(x1, list(map(spList[j],(x1))), colors[j])
  plt.subplot(base + 2*i + 1)
  #print(len(spList))
  for k in range(len(nVals)):
    iSpline = spline(funcList[i], makeEquiDistRange(nVals[k]))
    spError = composedError(funcList[i], iSpline)
    plt.plot(x1, list(map(spError,(x1))), colors[k])
plt.show()

### Part 2 of project begins here.
nVals = [5,10,15,20,25,30]

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

print(makePrettyTable(cbList, integral1, nVals, (0,2)))
print()
print(makePrettyTable(cbList,integral2,nVals,(0,1)))
print()
print()
print(makePrettyTable([gaussianQuad],integral1,gnVals,(0,2),names=["Gaussian"]))
print()
print(makePrettyTable([gaussianQuad], integral2,gnVals,(0,1),names=["Gaussian"]))
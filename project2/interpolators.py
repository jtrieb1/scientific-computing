import math
from operator import mul
from functools import reduce

def makeEquiDistRange(n):
  return [2*(i/n) - 1 for i in range(n+1)]

def makeCosRange(n):
  return [math.cos((2*i - 1)*math.pi / (2*n)) for i in range(1,n+1)]

def polynomialInterpolate(method, nodes):
  evalColumn = [method(nodes[i]) for i in range(len(nodes))]
  c = [[evalColumn[i]] + [0]*(len(nodes) - 1) for i in range(len(nodes))]
  for j in range(1, len(nodes)):
    for i in range(len(nodes)-j):
      c[i][j] = (c[i+1][j-1]-c[i][j-1])/(nodes[i+j]-nodes[i])
  
  def interpolated(x):
    total = 0
    for i in range(len(nodes)):
      total += c[0][i] * reduce(mul, [(x - nodes[j]) for j in range(i)], 1)
    return total
  
  return interpolated


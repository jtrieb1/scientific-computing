"""An implementation of Shepard's method of interpolation in higher dimensions."""
from operator import mul
from functools import reduce
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def dataParser(args):
  """Turns input into something the program can work with."""
  vals = args.split(";")
  vals = [val.strip() for val in vals]
  vals = [[float(x) for x in val[1:-1].split(",")] for val in vals]
  return vals

def squareDistance(value1, value2):
  """Standard non-negative evaluation function for Shepard's method."""
  return abs((value1[0]-value2[0])**2 + (value1[1]-value2[1])**2)
  
def shepard(vals,x):
  def v(i):
    funcList = [squareDistance(x,vals[j]) for j in range(len(vals)-1) if j != i]
    return reduce(mul, funcList)
  vsum = sum([v(i) for i in range(len(vals)-1)])
  def w(i):
    return v(i)/vsum
  return sum([vals[i][2]*w(i) for i in range(len(vals)-1)])

def main():
  data = input("Input your initial data points for interpolation in the format (x1,y1,c1);(x2,y2,c2);...:\n")
  values = dataParser(data)
  testPoints = input("Input your desired test values in the format (x1,y1);(x2,y2);...:\n")
  tests = dataParser(testPoints)
  print("Calculating interpolation...")
  output = [(point,shepard(values,point)) for point in tests]
  for line in output:
    print(str(line[0]) + ":  " + str(line[1]))
  
  # Now we generate a visual representation of what we've done.
  def contourFunc(x,y):
    return shepard(values,[x,y])
  x = np.linspace(min([val[0] for val in values])-4,max([val[0] for val in values])+4,50)
  y = np.linspace(min([val[1] for val in values])-4,max([val[1] for val in values])+4,50)
  X,Y = np.meshgrid(x,y)
  Z = contourFunc(X,Y)
  fig = plt.figure()
  ax = plt.axes(projection="3d")
  ax.plot_surface(X,Y,Z, rstride=1, cstride=1, cmap="viridis", edgecolor="none")
  ax.set_title("Interpolation map")
  plt.show()  
  

main()
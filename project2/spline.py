def zCalc(method, nodes):
  """This method solves the tridiagonal matrix associated to a given 
  function to yield the necessary z values for building the spline."""
  n = len(nodes) - 1
  yVals = list(map(method, nodes))
  h = [0]*(n)
  b = [0]*(n)
  for i in range(n):
    h[i] = nodes[i+1] - nodes[i]
    b[i] = 6*(yVals[i+1]-yVals[i])/h[i]
  u = [0]*(n)
  v = [0]*(n)
  u[1] = 2*(h[0]+h[1])
  v[1] = b[1]-b[0]
  for j in range(2,n):
    u[j] = 2*(h[j]+h[j-1]) - ((h[j-1]**2)/u[j-1])
    v[j] = b[j]-b[j-1]-(h[j-1]*v[j-1]/u[j-1])
  z = [0]*(n+1)
  z[-1] = 0
  for k in range(n-1,0,-1):
    z[k] = (v[k]-(h[k]*z[k+1]))/u[k]
  z[0] = 0
  return z, h, yVals

def spline(method, nodes):
  """This method actually builds the spline and returns a callable version of it,
  which will associate the correct values of x to their specific sub-function."""
  z, h, y = zCalc(method,nodes)
  A = [0]*len(nodes)
  B = [0]*len(nodes)
  C = [0]*len(nodes)
  for i in range(len(h)):
    A[i] = (1/(6*h[i]))*(z[i+1]-z[i])
    B[i] = z[i]/2
    C[i] = (-h[i]/6)*z[i+1] - (h[i]/3)*z[i] + (1/h[i])*(y[i+1]-y[i])

  def evaluate(x):
    for i in range(len(nodes)):
      if x > nodes[i] and x <= nodes[i+1]:
        return y[i] + (x - nodes[i])*(C[i] + (x-nodes[i])*(B[i] + (x-nodes[i])*A[i]))
  return evaluate



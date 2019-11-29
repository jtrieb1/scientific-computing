import math
def compositeTrapezoidal(method, n, interval):
  a = interval[0]
  b = interval[1]
  h = (b-a)/n
  xPartial = [a + i*h for i in range(1, n)]
  yPartial = list(map(method, xPartial))
  return (h/2)*(method(a) + 2*sum(yPartial) + method(b))

def compositeSimpson(method, n, interval):
  if n % 2 == 1:
    carry = True
  else:
    carry = False
  a = interval[0]
  b = interval[1]
  h = (b-a)/n
  x = [a + i*h for i in range(n+1)]
  presum = 0
  if carry:
    upperlimit = int((n-1)/2)
  else:
    upperlimit = int(n/2)
  innerTerm1 = [method(x[2*i - 2]) for i in range(2, upperlimit + 1)]
  innerTerm2 = [method(x[2*i - 1]) for i in range(1, upperlimit + 1)]
  presum = (h/3)*(method(a) + 2*sum(innerTerm1) + 4*sum(innerTerm2) + method(b))
  if not carry:
    return presum
  else:
    a = x[-2]
    finalTerm = ((b-a)/6)*(method(a) + 4*method((a+b)/2) + method(b))
    return presum + finalTerm

def legendreRootApprox(n,k):
  """This is an asymptotic approximation for the roots of the nth degree Legendre polynomial given by Tricomi. 
     I ended up not needing it, but how cool is it that this exists?"""
  return (1 - (1/(8*(n**2))) + (1/(8*(n**3))))*math.cos(math.pi*(4*k - 1)/(4*n + 2))

def gaussianQuad(method, n, interval):
  rootsDict = {
    2:  [0, (1/5)*math.sqrt(15)],
    4:  [0, (1/21)*math.sqrt(245-14*math.sqrt(70)), (1/21)*math.sqrt(245+14*math.sqrt(70))],
    8:  [0, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626],
    12: [0, 0.230458315955135, 0.448492751036447, 0.642349339440340, 0.801578090733310, 0.917598399222978, 0.984183054718588],
    16: [0, 0.178484181495848, 0.351231763453876, 0.512690537086477, 0.657671159216691, 0.781514003896801, 0.880239153726986, 0.950675521768768, 0.990575475314417],
    32: [0, 0.093631065854733, 0.186439298827992, 0.277609097152497, 0.366339257748073, 0.451850017272451, 0.533389904786348, 0.610242345836379, 0.681731959969743, 0.747230496449562, 0.806162356274167, 0.858009652676504, 0.902316767743434, 0.938694372611168, 0.966822909689993, 0.986455726230642, 0.997424694246455]
  }
  weightsDict = {
    2:  [(8/9),(5/9)],
    4:  [128/225, (1/900)*(322+13*math.sqrt(70)), (1/900)*(322-13*math.sqrt(70))],
    8:  [0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574],
    12: [0.232551553230874, 0.226283180262897, 0.207816047536889, 0.178145980761946, 0.138873510219787, 0.092121499837728, 0.040484004765316],
    16: [0.179446470356207, 0.176562705366993, 0.168004102156450, 0.154045761076810, 0.135136368468525, 0.111883847193404, 0.085036148317179, 0.055459529373987, 0.024148302868548],
    32: [0.093768446160210, 0.093356426065596, 0.092123986643317, 0.090081958660639, 0.087248287618844, 0.083647876067039, 0.079312364794887, 0.074279854843954, 0.068594572818657, 0.062306482530317, 0.055470846631664, 0.048147742818712, 0.040401541331670, 0.032300358632329, 0.023915548101749, 0.015321701512935, 0.006606227847587]
  }
  assert n in [2,4,8,12,16,32]
  a = interval[0]
  b = interval[1]
  u = (a+b)/2
  S = weightsDict[n][0]*method(u)
  for i in range(1,int(n/2)+1):
    u = ((b-a)*rootsDict[n][i] + a + b)/2
    v = ((a-b)*rootsDict[n][i] + a + b)/2
    S += weightsDict[n][i]*(method(u) + method(v))
  S = (b-a) * S / 2
  return S
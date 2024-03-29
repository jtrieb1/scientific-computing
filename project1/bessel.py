factorialList = [1]

def factorial(n):
    """ This is just a relatively fast, dynamic way of calculating factorials.
     Since I expect to call this function fairly often, saving the
     intermediate values cuts down computation time slightly at the
     expense of needing more memory to hold the array. """
    cur = len(factorialList)
    if cur <= n :
        for i in range(cur, n+1):
            factorialList.append(factorialList[i-1] * i)
    return factorialList[n]

def truncatedBessel(x,M):
    """ Calculates the Bessel function of order zero at the point x,
     using the truncated infinite series to obtain a reasonable
     approximation. """
    value = 0
    for p in range(M+1):
        innerCoeff = ((-1)**p)/(factorial(p)**2)
        variable = (x/2)**(2*p)
        value += innerCoeff * variable
    return value

def generalTruncatedBessel(x,n,M):
    """ Calculates the Bessel function of order n at the point x, using
        the infinite series truncated at M terms. """
    value = 0
    for p in range(M+1):
        innerCoeff = ((-1)**p)/(factorial(p)*factorial(n+p))
        variable = (x/2)**(n + 2*p)
        value += innerCoeff * variable
    return value

def dynamicBessel(x,N):
    """ Calculates the Bessel function of order zero at the point x,
     using the recursion formula and taking the appropriate corrective
     factor. Requires slightly more terms than truncatedBessel to achieve the
     same relative error bound. """
    J = [0]*(N+1)
    # It's important to pick N large enough so that the expression
    # $$sum_(p=0)^(infinity) [(((-1)**p)/(p!(n+p)!))(x/2)**(n+2p)]$$
    # becomes negligible for n >= N.
    J[N-1] = 1
    for i in range(1,N)[::-1]:
        intermediate = (2*(i)/x)*J[i] - J[i+1]
        J[i-1] = intermediate
    lambdaArray = [J[0]] + [2*x for x in J[2::2]]
    correctiveFactor = sum(lambdaArray)
    return J[0] / correctiveFactor

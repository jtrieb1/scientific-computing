def bisectionMethod(function, interval, maxIterations, epsilon, delta):
    """ Calculates and returns a root of the given function (if one exists)
        within the specified interval. The algorithm continues until either:
        1. maxIterations has been reached.
        2. The function applied to the current guess has modulus < epsilon.
        3. The associated error bound is less than delta.
    """
    a = float(interval[0])
    b = float(interval[1])
    u = function(a)
    v = function(b)
    sign = lambda x: 1 if x >= 0 else -1
    assert(sign(u) != sign(v)), "Bisection method inapplicable on this interval. Find an interval whose endpoints yield different signs."
    error_bound = b-a
    for k in range(1,maxIterations+1):
        # Bisect our interval
        error_bound = error_bound / 2.0
        # Get value at the bisection
        c = a + error_bound
        w = function(c)
        # If we are within our stated acceptable error, we are done
        if abs(error_bound) < delta or abs(w) < epsilon:
            return (c, w, k)
        # Otherwise, check if the root is between a and c or c and b
        if sign(w) != sign(u):
            # Root is between a and c, so check interval [a,c]
            b = c
            v = w
        else:
            # Root is between c and b, so check interval [c,b]
            a = c
            u = w
    # This code should be unreachable:
    return (0,0,-1)

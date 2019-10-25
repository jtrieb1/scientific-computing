def secantMethod(function, interval, maxIterations, epsilon, delta):
    """ Calculates a root of the given function within the given interval,
        using the secant method. """
    a = interval[0]
    b = interval[1]
    fa = function(a)
    fb = function(b)
    for k in range(2, maxIterations+1):
        if abs(fa) > abs(fb):
            # If this is the case, we need to swap a with b and fa with fb.
            tmp = a
            a = b
            b = tmp
            tmp = fa
            fa = fb
            fb = tmp
        s = (b - a)/(fb - fa)
        b = a
        fb = fa
        a = a - fa * s
        fa = function(a)
        if abs(fa) < epsilon or abs(b-a) < delta:
            return (a, fa, k)
    # This should be unreachable:
    return (0,0,-1)

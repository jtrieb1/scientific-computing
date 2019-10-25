def getMaxLength(matrix):
    return max([max([len(str(x)) for x in row]) for row in matrix])

def prettyTable(matrix, headers, xlabel, ylabel, sep="|", padding=1, hRule=True):
    assert(len(headers) + 1 >= max([len(row) for row in matrix])), "Not enough headers for data!"
    maxSpaceRequiredPerCell = getMaxLength(matrix) + padding
    yAxis = "{0:>{space}}".format(ylabel, space=maxSpaceRequiredPerCell) + sep
    xAxis = sep.join(["{0:^{space}}".format(xlabel, space=maxSpaceRequiredPerCell)] + [" "*maxSpaceRequiredPerCell]*len(headers)) + "\n"
    formattedHeaders = yAxis + sep.join(["{0:^{space}}".format(x, space=maxSpaceRequiredPerCell) for x in headers]) + "\n"
    formattedHeaders += xAxis
    horizontalLine = "{0:=^{space}}".format("",space=maxSpaceRequiredPerCell+1)*(len(headers)+1) + "\n"
    body = "\n".join([sep.join(["{0:^{space}}".format(x, space=maxSpaceRequiredPerCell) for x in row]) for row in matrix])
    return formattedHeaders + (horizontalLine if hRule else "") + body

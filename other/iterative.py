def permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    res = []
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return res
    indices = list(range(n))
    cycles = list(range(n, n-r, -1))
#    yield tuple(pool[i] for i in indices[:r])
    res.append(tuple(pool[i] for i in indices[:r]))
    #print(tuple(pool))
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
        #        yield tuple(pool[i] for i in indices[:r])
                res.append(tuple(pool[i] for i in indices[:r]))
                break
        else:
            return res

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    res = []
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return res
    indices = list(range(r))
    res.append(tuple(pool[i] for i in indices))
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return res
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        res.append(tuple(pool[i] for i in indices))

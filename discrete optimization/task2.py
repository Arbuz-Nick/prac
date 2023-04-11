from solver import ListSolver
from math import sqrt

def g(C, x, weights):
    return C - sum([weight * x_i for weight, x_i in zip(weights, x)])


def new_p(l, p, w):
    res = []
    for i in range(len(p)):
        tmp = p[i]
        for j in range(len(l)):
            tmp -= l[j] * w[i][j]
        res.append(tmp)
    return res


def D(x, l, p, C):
    res = 0
    # print(p)
    for i in range(len(l)):
        res += l[i] * C[i]
    #     print(i, res)
    res += sum([x[i] * p[i] for i in range(len(x))])
    return res

def norm(l):
    return sqrt(sum([i**2 for i in l]))
n, m = map(int, input().split())

C = list(map(int, input().split()))
#for i in range(m):
#    C.append([i])

p = list(map(int, input().split()))
#for i in range(n):
#    p.append(l[i])

w = []
for i in range(n):
    w.append(list(map(int, input().split())))

l = []
for i in range(m - 1):
    l.append(1)
prev_D = 0
new_D = 0
for k in range(100000):
    # print()
    new_Ps = new_p(l, p, [j[:-1] for j in w])
    new_Ws = [i[-1] for i in w]
    # print(new_Ws)
    # print([p_ for p_ in new_Ps if p_ > 0])
    sol = ListSolver(n, C[-1])
    max_S, indexes = sol.solve(new_Ws, new_Ps)
    #print(max_S, indexes)
    x = []
    for i in range(n):
        if i in indexes:
            x.append(1)
        else:
            x.append(0)
    new_D = D(x, l, new_Ps, C[:-1])
    delta_D = []
    for i in range(m - 1):
        delta_D.append(g(C[i], x, [j[i] for j in w]))
    norm_D = norm(delta_D)
    #if (k % 10000 == 0):
    #print(l, delta_D, norm_D, k, new_D)
    for i in range(m - 1):
        if norm_D > 0:
            l[i] = max(l[i] - delta_D[i] / norm_D / (k + 1), 0)
    if abs(new_D - prev_D) < 1e-6:
        break
    prev_D = new_D
print(new_D)
for i in l:
    print(i)
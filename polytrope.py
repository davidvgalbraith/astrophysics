import math 

#Polytrope approximation problem, by David Galbraith

def approximate():
    rSun = 70000000000.0
    G = 7 * (10**-8)
    K = 6 * 10**13
    Q = []
    P = []
    P.append(150.0);
    Q.append(10000000.0);
    iterations = 100000
    deltar = rSun/iterations
    r = .007 * rSun
    for i in range(1, iterations):
        r += deltar
        Q.append(Q[i-1] - deltar * P[i-1] * r * r * G * 4 * math.pi * 5 /(6 * K))
        dpdr = Q[i] * P[i-1] ** .8 / (r * r)
        P.append(P[i-1] + dpdr * deltar)
        print P[i]
    

approximate()

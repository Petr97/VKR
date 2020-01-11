"""
Функция, вычисляющая С(m,n)*a^m mod p

Вход  : n,m,a
Выход : С(m,n)*a^m mod p
"""
def mod_binomial_x_a(n,m,a):
    K = a.base_ring()
    if n<m:
        print "n<m"
        return 0
    if m == 0 or n == 0:
        return K(0)
    s1,s2 = K(1),K(1)
    for i in range(n-m+1,n+1):
        s1 = (s1*i)
    for i in range(2,m+1):
        s2 = (s2*i)
    return s1*(a^m)/s2


"""
Функция, возращающая произвольную точку кривой.

Вход  : кривая C
Выход : случайная точка на кривой С
"""
def random_point(C):
    out = []
    f = C.hyperelliptic_polynomials()[0]
    while(True):
        X = (C.base_ring()).random_element()
        if (f(X).is_square()):
            Y=sqrt(f(X))
            out=[X,Y]
            return C(out)


"""
Функция, возводящая коэффициенты матрицы А в степень р.

Вход  : матрица А, р
Выход : А^(p)
"""
def pow_matrix(A,p):
    B = matrix(A)
    for i in range(0,2):
        row = [B[i][0]^p,B[i][1]^p]
        B.set_row(i,row)
    return B


"""
Функция, реализующая метод Картье-Манина

Вход  : f(x) - многочлен кривой, 
        р - характеристика поля, n - степень расширения
Выход : (s1,s2) mod p
"""
def manin(f,p,n):
    if len(f.coefficients())==2 and f[0]!=0 and f[5]!=0:
        z = int((p-1)/2)
        k1,k2 = 3*z//5 if z%5==0 else 0, (3*z+1)//5 if (3*z+1)%5==0 else 0
        k3,k4 = (z-1)//5 if (z-1)%5==0 else 0, z//5 if z%5==0 else 0
        A = matrix([[mod_binomial_x_a(z,k1,f[0]),mod_binomial_x_a(z,k2,f[0])],\
                    [mod_binomial_x_a(z,k3,f[0]),mod_binomial_x_a(z,k4,f[0])]])
    else:
        fp = f^((p-1)/2)
        A = matrix([[fp[p-1],fp[p-2]],[fp[2*p-1],fp[2*p-2]]])
    Af = matrix(A)
    for i in range(1,n):
        Af = Af*pow_matrix(A,p^(i))
    s1 = Af[0][0]+Af[1][1]
    s2 = Af[0][0]*Af[1][1]-Af[1][0]*Af[0][1]
    return (s1,s2)


"""
Функция, восстанавливающая многочлен х(t) с помощью
множества случайных дивизоров (над простым полем)

Вход  : hi_mod - (s1,s2) mod p, 
        р - характеристика поля, 
Выход : (s1,s2)
"""
def recovery_hi(hi_mod,p,f):
    s1,s2 = 0,int(hi_mod[1])
    bound_s1 = int(4*sqrt(float(p)))
    if int(hi_mod[0]) > bound_s1 or int(hi_mod[0]) < -bound_s1:
        s1 = -int(-hi_mod[0])
    else:
        s1 = int(hi_mod[0])
    bound_s2_1 = -int(6*p)
    bound_s2_2 = int(6*p)
    candidates_s2 = [s2]
    print "{} <= s2 <= {}".format(bound_s2_1,bound_s2_2)
    i = p
    while (s2-i >= bound_s2_1):
        candidates_s2.append(s2-i)
        i+=p
    i = p
    while (s2+i <= bound_s2_2):
        candidates_s2.append(s2+i)
        i+=p
    C = HyperellipticCurve(f,0)
    J = C.jacobian()
    while len(candidates_s2)!=1:
        D = J(random_point(C))
        for j in candidates_s2:
            ord_J = 1-s1+j-s1*p+p^2
            if not (D*int(ord_J)).is_zero():
                candidates_s2.remove(j)
    s2 = candidates_s2[0]
    return (s1,s2)


"""
Тест
"""
from time import time

p,n = 111999991,1
K = GF(p)                                               
R.<x> = PolynomialRing(K)
f = x^5+123*x^4+54*x^3+6423*x^2+625*x+467
f = x^5+46714179
tm1 = time()
hi_mod = manin(f,p,n)
print "коэффициенты x(t) mod p: {}".format(hi_mod)
print""
hi = recovery_hi(hi_mod,p,f)

print "Затраченное время: {}  c.".format(time()-tm1)
print ""
s1,s2 = hi
ord_J = 1-s1+s2-s1*p+p^2
print "x(t) = t^4 - ({})t^3 + ({})t^2 - ({})t + ({})".format(s1,s2,s1*p,p^2)
print "#J = {}".format(ord_J)



"""
проверка "в лоб" (через подсчет точек кривой)
"""

C=HyperellipticCurve(f,0)
CFq,CFq2 = C.count_points(2)
s1 = p-CFq+1
s2 = (CFq2-p^2 +s1^2-1)/2
ord_J = 1-s1+s2-s1*p+p^2
print "x(t) = t^4 - ({})t^3 + ({})t^2 - ({})t + ({})".format(s1,s2,s1*p,p^2)
print "#J = {}".format(ord_J)

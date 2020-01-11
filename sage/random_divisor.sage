"""
Функция, возращающая случайную точку кривой.

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
Функция, осуществляющая бинарный поиск элемента [a,r] из массива А = [[a0,r0]...]
(поиск осуществляется по элементу a). Если он не находит такого a,
он вставляет в массив элемент [a,r] на место (не нарушая сортировки)

Вход  : элемент r,a и массив A
Выход : позиция элемента [a,r] или -1
"""
def bin_search_insert(r,a,A):
    i,b_left,b_right = len(A)//2,0,len(A)-1
    while True:
        if b_left > b_right:
            A.insert(i+1,[a,r])
            return -1
        elif a>A[i][0]:
            b_left = i+1
            i = (b_right+b_left)//2
        elif a<A[i][0]:
            b_right = i-1
            i = (b_right+b_left)//2
        else:
            return i



"""
Функция, находящая элементы массива А = [[a0,b0,c0]...] для которых
ai == aj и bi == bj и ci != cj. Функция складывает такие элементы в
массив out (в одном экземпляре).

Вход  : массив А
Выход : массив out
"""
def search_same(A):
    out = []
    for i in range(0,len(A)-1):
        for j in range(i+1,len(A)):
            if A[i][0][0]==A[j][0][0] and A[i][0][1]==A[j][0][1] and A[i][2]!=A[j][2]:
                out.append(abs(A[i][1]-A[j][1]))
    if len(out)==0:
        return (0,0)
    return out



"""
Функция, котрая должна выполнятся параллельно (см 3.2)

Вход  : i - номер потока, w - ширина границы Хассе-Вейля,
        h - хеш-функция, D - дивизор, E = с*D ,
        li - массив шагов, Di - массив li*D,
        r - случайное целое r < w, b - случайный бит
        доп.параметр ncpus = M
Выход : столкнувшиеся элементы ( [R1,R2],r,b ),
        где R1,R2 -- координаты столкнувшихся дивизоров 
        в представлении Мамфорда
"""
@parallel(ncpus=16)
def for_treading(i,w,h,D,E,li,Di,r,b):
    R = r*D + b*E
    Rr = [[R,r]]
    j = 1
    while True:
        r += li[h(R)]
        R += Di[h(R)]
        j+=1
        if bin_search_insert(r,R,Rr)!=-1:
            print "количиство шагов до столкновения в {}-ом потоке {}:".format(i,j)
            return ([R[0],R[1]],r,b,j)



"""
Функция, осуществляющая распаралеливание

Вход  : w - ширина границы Хассе-Вейля,
        h - хеш-функция, D - дивизор, E = с*D ,
        li - массив шагов, Di - массив li*D,
        M - количество потоков, 
        global_Ri - массив столкнувшихся элементов
Выход : --
"""
def treading(w,h,D,E,li,Di,M,global_Ri):
    tri = [(i,w,h,D,E,li,Di,randint(0,w),randint(0,1)) for i in range(0,M)]
    res = for_treading(tri)
    for i in res:
        a = i[1]
        global_Ri.append(a)


"""
Функция, осуществляющая поиск порядка дивизора D

Вход  : C -  кривая, D - дивизор
Выход : порядок D или -1 при неудаче
"""
def order_divisor(C,D,Mg,M):
    q = C.base_ring().order()
    w = 2*int(4*(q+1)*sqrt(float(q)))
    c = q^2+6*q+1
    global_Ri = []
    h = lambda x:(int(x[0][0]))%Mg
    bMsq = M*int(sqrt(float(w)))
    li = [2*randint(int(sqrt(float(bMsq))),bMsq) for i in range(0,Mg)]
    Di = [D*li[i] for i in range(0,Mg) ]
    E = c*D
    treading(w,h,D,E,li,Di,M,global_Ri)
    avg_j = 0
    for i in global_Ri:
        avg_j+=i[-1]
    print "среднее количество шагов в потоках: {}".format(float(avg_j)/M)
    r1r0 = search_same(global_Ri)
    if r1r0==(0,0):
        return []
    n = [c-r1r0[i] for i in range(0,len(r1r0))]
    out = []
    for i in n:
        if (i*D).is_zero() and i!=0:
            out.append(abs(i))
    if len(out) == 0:
        return -1
    n = gcd(out)
    fd = list(factor(n))
    for i in range(0,len(fd)):
        for j in range(0,fd[i][1]):
            if not ((n//(fd[i][0]^j))*D).is_zero():
                n = n//(fd[i][0]^(j-1))
                break
    print "sqrt(w)/M = {}".format(sqrt(float(w)))
    return n



"""
Функция, осуществляющая поиск претендентов порядка якобиана
из порядка произвольного дивизора D

Вход  : C -  кривая, out - порядок
Выход : массив претендентов на порядок якобиана

"""

def recovery_N(out,C):
    q = C.base_ring().order()
    w = 2*int(4*(q+1)*sqrt(float(q)))
    n = out
    
    bound1 = int( (sqrt(float(q))-1 )**4+0.5 )
    bound2 = int( (sqrt(float(q))+1 )**4-0.5 )
    j = 1
    while not bound1<=n<=bound2:
        n+=out
        j+=1
    print "количество умножений, что бы достичь границы Хассе-Вейля: {}".format(j)
    candidates = []
    
    variants_N = 1 if n>=w else round(float(w+1)/n+0.5)
    for i in range(0,variants_N+1):
        candidates.append((i+1)*out)
    
    return candidates


"""
Тест
"""
from time import time

p,n = next_prime(19001),1
K = GF(p)                                               
R.<x> = PolynomialRing(K)
f = x^5+123*x^4+54*x^3+6423*x^2+625*x+467


Mg,M = 10,16                                # отличительная характеристика и число потоков
C = HyperellipticCurve(f,0)
J = C.jacobian()
D = J(random_point(C))
print D

tm1 = time()
out = order_divisor(C,D,Mg,M)
if out!=-1:
    print "порядок дивизора D = {}: {}".format(D,out)
    N = recovery_N(out,C)
    print "кондитаты в #J: {}".format(N)
else:
    print "неудача!"

print "Затраченное время: {}  c.".format(time()-tm1)



"""
проверка "в лоб" (через подсчет точек кривой)
"""

CFq,CFq2 = C.count_points(2)
s1 = p-CFq+1
s2 = (CFq2-p^2 +s1^2-1)/2
ord_J = 1-s1+s2-s1*p+p^2
print "проверка с помощью посчета точек: "
print "x(t) = t^4 - ({})t^3 + ({})t^2 - ({})t + ({})".format(s1,s2,s1*p,p^2)
print "#J = {}".format(ord_J)



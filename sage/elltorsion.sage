"""
Блок алгоритма интерполяции (подробнее см cf Algo 10.3,
page 281 of Modern Computer Algebra)

главная функция interpolate(xi,yi)
Ввод  : xi - аргументы функции 
        yi=f(xi) - значения функции
Вывод : f - результат интерполяции
"""
############################################################################################################################
def build_heap(xi):
    k = 0
    while ((1<<k)<len(xi))!=0:
        k+=1
    m = 2*(1<<k)
    heap = [1 for i in range(0,m)]
    for i in range(0,len(xi)):
        heap[(m>>1)+i] = x-xi[i]

    m >>= 2
    while (m > 0):
        for i in range(0,m):
            heap[m+i] = heap[2*(m+i)]*heap[2*(m+i)+1]
        m >>= 1
    return heap
    
def multieval(f,heap,n):
    m = len(heap)
    v,w,yi = [f%heap[1]],[],[]
    i = 2
    while i<m:
        for j in range(0,i//2):
            w.append(v[j]%heap[i+2*j])
            w.append(v[j]%heap[i+2*j+1])
        v = w
        w = []
        i<<=1
    for i in range(0,n):
        if v[i].degree()!=0:
            print "err in multieval"
            1/0
        yi.append(v[i])
    return yi

def interpolate(xi,yi):
    sub_prod = build_heap(xi)
    if len(xi)!=len(yi):
        print " len(xi) != len(yi) "
        1/0
    pol = diff(sub_prod[1])
    si = multieval(pol,sub_prod,len(xi))
    for i in range(0,len(yi)):
        si[i] = yi[i]//si[i]
    m = len(sub_prod) >> 1
    v = [0 for i in range(0,m)]
    for i in range(0,len(yi)):
        v[i] = si[i]
    m>>=1
    while m>0:
        for i in range(0,m):
            v[i] = v[2*i]*sub_prod[2*m+2*i+1]+v[2*i+1]*sub_prod[2*m+2*i]
        m>>=1
    return v[0]
############################################################################################################################

"""
Блок алгоритма поиска многочленов деления Кантора
(подробнее см *годри1* )

главная функция Cantor_Poly(f,n)
Ввод  : f - многочлен кривой 
        n - размер кручения
Вывод : (d0,d1,d2,e0,e1,e2) - многочлены деления кантора
"""
############################################################################################################################

temp_value_1 = [None,[0,0]]                    
def __psi(f,n):
    if temp_value_1[0]!=f:                     
        temp_value_1[0] = f
        temp_value_1[1] = [0,0]
    if n == 0 or n == 1:
        return 0
    if len(temp_value_1[1])>n and temp_value_1[1][n]!=0:
        return temp_value_1[1][n]
    if n == 2:
        psi = 1
    elif n == 3:
        psi = 4*f
    elif n == 4:
        p1 = -f[3]-4*f[4]*x-10*x**2
        p2 = (-3*f[3]*x**2-5*x**4-4*f[4]*x**3-2*f[2]*x-f[1])*(3*f[3]*x+f[2]+6*f[4]*x**2+10*x**3)
        p3 = (-3*f[3]*x**2-5*x**4-4*f[4]*x**3-2*f[2]*x-f[1])**3
        psi = -2*(p3-4*f*(p2-2*p1*f))
    elif n == 5:
        fp = -diff(f)
        p1 = f[4]+5*x
        p2 = -f[3]-4*f[4]*x-10*x**2
        p3 = (3*f[3]*x+f[2]+6*f[4]*x**2+10*x**3)
        p4 = p3**2
        psi = -4*f*(-5*fp**4+8*f*(3*p3*fp**2+2*f*(-p4-2*fp*p2+4*p1*f)))
    elif n == 6:
        p1 = -diff(f)
        p2 = 3*f[3]*x+f[2]+6*f[4]*x**2+10*x**3
        p3 = -(10*x**2+4*f[4]*x+f[3])
        p4 = f[4]+5*x
        psi = -4*(8*p3*f**2-4*p1*f*p2+p1**3)*(128*f**4+64*p2*p3*f**3+64*f**3*p1*p4-48*p3*(f*p1)**2- \
                48*f**2*p1*p2**2+40*p1**3*p2*f-7*p1**5)-(-32*f**2*p3*p1-5*p1**4+24*f*p1**2*p2-16*(p2*f)**2+64*p4*f**3)**2
    elif n == 7:
        p1 = -3*f[3]*x**2-5*x**4-4*f[4]*x**3-2*f[2]*x-f[1]
        p2 = 3*f[3]*x+f[2]+6*f[4]*x**2+10*x**3
        p3 = -f[3]-4*f[4]*x-10*x**2
        p4 = f[4]+5*x
        psi = -8*f*(256*f**4*p4*p2-256*f**4*p1+128*f**4*p3**2-192*f**3*p4*p1**2-384*p2*p3*f**3*p1-64*f**3*p2**3+ \
                160*p3*f**2*p1**3+240*p2**2*f**2*p1**2-140*p2*f*p1**4+21*p1**6)* \
                (64*f**3*p4-5*p1**4-32*f**2*p1*p3+24*p2*f*p1**2-16*f**2*p2**2)-16*f* \
                (64*f**3*p4*p1-7*p1**5+40*f*p1**3*p2-48*f**2*p1**2*p3-48*f**2*p2**2*p1+128*f**4+64*f**3*p2*p3)**2
    else:
        r = n//2
        s = n-r
        if r==s:
            r-=1
            s+=1
        if r==s-1:
            r-=1
            s+=1
        
        m11,m12,m13 = __psi(f,s)*__psi(f,r-2),__psi(f,s+1)*__psi(f,r-1),__psi(f,s+2)*__psi(f,r)
        m21,m22,m23 = __psi(f,s-1)*__psi(f,r-1),__psi(f,s)*__psi(f,r),__psi(f,s+1)*__psi(f,r+1)
        m31,m32,m33 = __psi(f,s-2)*__psi(f,r),__psi(f,s-1)*__psi(f,r+1),__psi(f,s)*__psi(f,r+2)
        
        det_m = m11*(m22*m33-m23*m32)+m12*(-m21*m33+m23*m31)+m13*(m21*m32-m22*m31)
        d = __psi(f,s)*__psi(f,r)*__psi(f,s-r)
        psi = -det_m // d                                                                 
        
    while len(temp_value_1[1])<=n:
        temp_value_1[1].append(0)
    temp_value_1[1][n] = psi
    return psi

temp_value_2 = [None,[1,1]]                    
def __gamma(f,n):
    if temp_value_2[0]!=f:                     
        temp_value_2[0] = f
        temp_value_2[1] = [1,1]
    if n == 0 or n == 1:
        return 1
    if len(temp_value_2[1])>n and temp_value_2[1][n]!=0:
        return temp_value_2[1][n]
    p1 = -3*f[3]*x**2-5*x**4-4*f[4]*x**3-2*f[2]*x-f[1]
    p2 = 3*f[3]*x+f[2]+6*f[4]*x**2+10*x**3
    p3 = -f[3]-4*f[4]*x-10*x**2
    p4 = f[4]+5*x
    
    if n == 2:
        gamma = 4*f
    elif n == 3:
        gamma = -64*p4*f**3+32*p1*p3*f**2+16*p2**2*f**2-24*p1**2*p2*f+5*p1**4
    elif n == 4:
        gamma = 1024*f**5+512*f**4*p1*p4+512*f**4*p2*p3-384*f**3*p1**2*p3-384*f**3*p1*p2**2+320*f**2*p1**3*p2-56*f*p1**5
    elif n == 5:
        gamma = -4096*f**6*p2**2+2048*f**5*p1**2*p2+16384*f**7*p4-256*f**4*p1**4-4096*f**6* \
        p3**3-14*p1**9-2560*f**4*p3**2*p1**3-384*f**2*p3*p1**6+512*f**4*p1*p2**4+512*f**3*p1**3*p2**3- \
        576*f**2*p1**5*p2**2+160*f*p1**7*p2-768*f**3*p1**5*p4+8192*f**6*p4**2*p1-4096*f**5*p3*p1**2*p4+ \
        10240*f**5*p3**2*p1*p2-6144*f**4*p3*p1**2*p2**2+3072*f**3*p3*p1**4*p2-4096*f**5*p1*p2**2*p4+4096*f**4*p1**3*p2*p4
    elif n == 6:
        gamma = 524288*f**9*p1-262144*f**9*p3**2-40960*f**5*p1**6-288*f*p1**11-327680*f**7*p1**3*p3- \
        393216*f**7*p1**2*p2**2+262144*f**8*p1**2*p4+262144*f**6*p1**4*p2-131072*f**7*p4**2*p1**3- \
        524288*f**9*p4**2*p3+16384*f**4*p4*p1**7-36864*f**5*p1**5*p3**2-98304*f**7*p1**2*p3**3- \
        5248*f**3*p1**8*p3+36864*f**4*p2**3*p1**5-40960*f**5*p2**4*p1**3+32768*f**6*p2**5*p1- \
        32768*f**7*p2**4*p3-17408*f**3*p2**2*p1**7+3712*f**2*p1**9*p2-131072*f**8*p2*p3**3+ \
        786432*f**8*p1*p2*p3-147456*f**5*p4*p1**5*p2+180224*f**6*p4*p1**4*p3+393216*f**6*p4*p1**3*p2**2- \
        786432*f**7*p4*p1**2*p2*p3-262144*f**7*p4*p2**3*p1+524288*f**8*p4**2*p1*p2+524288*f**8*p4*p3**2*p1+ \
        262144*f**8*p4*p3*p2**2+38912*f**4*p1**6*p3*p2-61440*f**5*p1**4*p3*p2**2+98304*f**6*p1**3*p3**2*p2- \
        32768*f**6*p1**2*p3*p2**3+196608*f**7*p1*p3**2*p2**2
    else:
        P = __psi(f,n+1)*__gamma(f,n-1)+__psi(f,n-1)*__psi(f,n+2)
        gamma = P//__psi(f,n)                                                            
        
    while len(temp_value_2[1])<=n:
        temp_value_2[1].append(0)
    temp_value_2[1][n] = gamma
    return gamma

temp_value_3 = [None,[-2,-2]]                    
def __alpha(f,n):
    if temp_value_3[0]!=f:                     
        temp_value_3[0] = f
        temp_value_3[1] = [-2,-2]
    if n == 0 or n == 1:
        return -2
    if len(temp_value_3[1])>n and temp_value_3[1][n]!=0:
        return temp_value_3[1][n]
    p1 = -3*f[3]*x**2-5*x**4-4*f[4]*x**3-2*f[2]*x-f[1]
    p2 = 3*f[3]*x+f[2]+6*f[4]*x**2+10*x**3
    p3 = -f[3]-4*f[4]*x-10*x**2
    p4 = f[4]+5*x
    if n == 2:
        alpha = 0
    elif n == 3:
        alpha = -2*p1
    elif n == 4:
        alpha = -8*f*p1
    elif n == 5:
        alpha = 64*p1*p3*f**2-40*p1**2*p2*f+9*p1**4-64*p4*f**3+16*p2**2*f**2
    elif n == 6:
        alpha = 1024*f**4*p1*p4-640*f**3*p1**2*p3-512*f**3*p1*p2**2+512*f**2*p1**3*p2-96*f*p1**5+1024*f**5+512*f**4*p2*p3
    else:
        P = __psi(f,n-1)*__alpha(f,n-1)+__psi(f,n)*__psi(f,n-3)
        alpha = P//__psi(f,n-2)
        
    while len(temp_value_3[1])<=n:
        temp_value_3[1].append(0)
    temp_value_3[1][n] = alpha
    return alpha

def __delta(f,n):
    d0 = -__psi(f, n-1)*__psi(f, n+1)
    d1 = -__psi(f, n-1)*__gamma(f, n)+__psi(f, n+1)*__alpha(f, n)
    d2 = -16*(f*__psi(f, n))**2
    return (d0,d1,d2)

def __epsilon(f,n):
    d0,d1,d2 = __delta(f,n)
    d0_p,d1_p,d2_p = __delta(f,n-1)
    d0_n,d1_n,d2_n = __delta(f,n+1)
       
    e2 = -__psi(f,n+1)*__psi(f, n-1)*(__psi(f, n)*d2)**2
    
    pp2 = (__psi(f, n-1))**2
    pn2 = (__psi(f, n+1))**2
    
    e0 = (-d2*pp2*d1_n+d2*pn2*d1_p+d1*pp2*d2_n-d1*pn2*d2_p)*d0
    e1 = d2**2*pp2*d0_n-d2**2*pn2*d0_p-d2*d0*pp2*d2_n+d2*d0*pn2*d2_p-d1*d2*pp2*d1_n+ \
    d1*d2*pn2*d1_p+d1**2*pp2*d2_n-d1**2*pn2*d2_p
    g = gcd(e1,gcd(e0,e2))
    while g.degree() >= 1:
        e0 //= g
        e1 //= g
        e2 //= g
        g = gcd(e1,gcd(e0,e2))
    return (e0,e1,e2)

def Cantor_Poly(f,n):
    d0,d1,d2 = __delta(f,n)
    d2 //= f
    d0 = x**2*d2+4*x*d1+16*d0*f
    d1 =-2*(2*d1 + x*d2)                        
    e0,e1,e2 = __epsilon(f,n)
    e0 = 4*f*e0+x*e1
    e1 = -e1
    e2 = -4*f*e2                                                                     
    return (d0,d1,d2,e0,e1,e2) 

############################################################################################################################

"""
Блок алгоритма поиска произведения вида [х-(q_i + r_j)],
где q_i и r_j - корни многочленов Q и R
(подробнее см Bostan, A., Flajolet, P., Salvy, B.,
Schost, 2006. Fast computation of special resultants.
Journal of Symbolic Computation 41,1–29.)

главная функция res_sum_roots(qQ,rR)
Ввод  : qQ, rR - унитарные многочлены
Вывод : искомое произведение [х-(q_i + r_j)]
"""
############################################################################################################################
def logTrunc(a,n):
    if (a.base_ring().characteristic()<n):
        print "err characteristic is to small log_trunc"
    if a[0] != 0:
        print "non-zero constant log_trunc"
    tm = diff(a)*((1+a).inverse_series_trunc(n))
    tm = a.parent()([tm[i]/(i+1) for i in range(0,n)])
    tm*=x
    return tm%(x**n)

def expTrunc(a,n):
    if (a.base_ring().characteristic()<n):
        print "err characteristic is to small log_trunc"
    if a[0] != 0:
        print "non-zero constant log_trunc"
    prec,r = 1,a.parent(1)
    while(True):
        prec<<=1
        tmp = 1+a-logTrunc(r-1,prec)
        r = r.multiplication_trunc(tmp,prec)
        if prec>=n:
            break
    return r%(x**n)

def res_sum_roots(qQ,rR):
    d = qQ.degree()*rR.degree()
    if d == 0:
        return 1
    Q,R = (qQ.monic()).reverse(),(rR.monic()).reverse()
    prec = d+3
    Qser = diff(Q).multiplication_trunc(Q.inverse_series_trunc(prec),prec)
    Rser = diff(R).multiplication_trunc(R.inverse_series_trunc(prec),prec)
    n = mod(2,Q.base_ring().characteristic())                                             
    i = 1
    Qt,Rt = Qser.coefficients(),Rser.coefficients()
    #print (Qt,Rt)
    while (i<prec):
        Qt[i] /= n
        Rt[i] /= n                                                       
        n*=(i+2)
        i+=1
    Qt.insert(0,-(qQ.degree()))
    Rt.insert(0,-(rR.degree()))
    Qser,Rser = qQ.parent()(Qt),qQ.parent()(Rt)
    
    Sser = Qser.multiplication_trunc(Rser,prec)
    n = mod(1,Q.base_ring().characteristic())
    i = 1
    St = Sser.coefficients()
    while (i<prec):
        St[i] = -St[i]*n
        n*=i
        i+=1
    St[0] = 0
    Sser = qQ.parent()(St)
    Sser = expTrunc(Sser,prec)
    return Sser.reverse(d)

############################################################################################################################

"""
Функция для поиска предпоследнего ненулевого подрезультанта

Ввод  : A, B - унитарные многочлены
        
Вывод : r1 - искомый подрезультант
"""
def subresultants(A, B):
    r1,r2 = A,B
    first = True
    while True:
        d1 = r1.degree()-r2.degree()
        y1 = r2.leading_coefficient()
        if first:
            B1 = (-1)^(d1+1)
            f1 = y1.parent()(-1)
            first = False
        else:
            B1 = -y0*f1^d1
            f1 = ((-y1)^(d1-1))/(f1^(d0-1))
        r = ((y1^(d1+1)*r1)%r2)/B1
        if r==0:
            return r1
        r1=r2
        r2=r
        d0=d1
        y0=y1
    
"""
Функция для разбиения многочлена на два многочлена
с четными и нечетными коэффициентами

Ввод  : многочлен f
        
Вывод : f_odd,f_even - многочлены с нечетными
и четными коэффициентами соответственно
"""
def split_odd_even(f):
    f_odd,f_even = 0,0
    for i in range(0,f.degree()//2):                         
        f_odd+=(f[2*i+1]*x**i)
        f_even+=(f[2*i]*x**i)
    if f.degree()%2!=0:
        f_odd+=f[f.degree()]*(x**(f.degree()//2))
    return (f_odd,f_even)



"""
Функция для вычисления многочленов G1(s,x) и G2(s,x)

Ввод  : d0,d1,d2 - многочлены деления Кантора
        s - текущая константа s
        
Вывод : многочлены G1 и G2
"""
def eval_G(d0,d1,d2,s):
    B = d0(x+s/2)
    Bodd,Beven = split_odd_even(B)
    Bodd,Beven = Bodd(-x+(s**2)/4),Beven(-x+(s**2)/4)
    dU0 = [Beven - (s/2)*Bodd, Bodd]
    B = d1(x+s/2)
    Bodd,Beven = split_odd_even(B)
    Bodd,Beven = Bodd(-x+(s**2)/4),Beven(-x+(s**2)/4)
    dU1 = [Beven - (s/2)*Bodd, Bodd]
    B = d2(x+s/2)
    Bodd,Beven = split_odd_even(B)
    Bodd,Beven = Bodd(-x+(s**2)/4),Beven(-x+(s**2)/4)
    dU2 = [Beven - (s/2)*Bodd, Bodd]

    
    G1 = dU0[1]*dU2[0] - dU0[0]*dU2[1]                           
    G2 = dU1[1]*dU2[0] - dU1[0]*dU2[1]

    return (G1,G2)



"""
Функция для вычисления паразитов

Ввод  : d2 - многочлен деления Кантора
        f - многочлен кривой
        si - множество констант
        
Вывод : parasites - паразиты для результанты
        subparasites - подпаразиты для подрезультант
"""
def eval_parasites(d2,f,si):                                                     
    if (d2%(f**3)!=0):
        print "computeParasites: d2 is not divisible by f^3"
    fac_d2 = d2//(f**3)
    fac_d2 = gcd(fac_d2,diff(fac_d2))
    
    double_d2 = (fac_d2(x/2)).monic()
    sum_roots_d2 = res_sum_roots(fac_d2,fac_d2)
    sum_roots_d2 //= double_d2
    sum_roots_d2 = gcd(sum_roots_d2,diff(sum_roots_d2))
    
    double_f = (f(x/2)).monic()
    sum_roots_f = res_sum_roots(f,f)
    sum_roots_f //= double_f
    sum_roots_f = gcd(sum_roots_f,diff(sum_roots_f))
    cross = res_sum_roots(fac_d2,f)
    ro = cross^6*double_d2*sum_roots_d2^4*double_f^3*sum_roots_f^9
    sro = cross^2*sum_roots_d2*sum_roots_f^4
    
    parasites = [0 for i in si]
    subparasites = [0 for i in si]
    for i in range(0,len(si)):
        parasites[i] = ro(si[i])
        subparasites[i] = sro(si[i])
    return parasites,subparasites



"""
Функция для вычисления множеств [результанта/паразит]
и множеств [подрезультанта0/подпаразит] и [подрезультанта1/подпаразит]

Ввод  : d0,d1,d2 - многочлены деления Кантора
        si - множество констант
        parasites - множество паразитов
        sub_parasites - множество подпаразитов
        
Вывод : res,subres0,subres1 - искомые множества
"""
def eval_res(d0,d1,d2,si,parasites,sub_parasites):
    res = [0 for i in si]
    subres0 = [0 for i in si]
    subres1 = [0 for i in si]
    for i in range(0,len(si)):
        G1,G2 = eval_G(d0,d1,d2,si[i])
        G1,G2 = G1.monic(),G2.monic()
        res[i] = (G1.resultant(G2))/parasites[i]                   
        subres = subresultants(G1,G2)
        subres0[i],subres1[i] = subres[0]/sub_parasites[i],subres[1]/sub_parasites[i]
        
    return (res,subres0,subres1)

"""
Функция для вычисления идеала кручения (недоделанная)
Ввод  : f - многочлен кривой
        l - размер скручивания
        
Вывод : res - ~R1(x)
        param - -S1/S2

недоделанно:
        G3 и НОД(G3,~R1)
        V0/V1
        V1
"""
def ell_torsin_ideal(f,l):
    d0,d1,d2,e0,e1,e2 = Cantor_Poly(f,l)
    n = ( (l**2*(7*l**2 - 33)) / 2 ) + 22 + 5
    si = [f.base_ring()(i+1) for i in range(0,n)]
    
    
    parasites,sub_parasites = eval_parasites(d2,f,si)
    si_res,si_sub_res0,si_sub_res1 = eval_res(d0,d1,d2,si,parasites,sub_parasites)
    res = interpolate(si,si_res).monic()
    sub_res0 = interpolate(si,si_sub_res0).monic()
    sub_res1 = interpolate(si,si_sub_res1).monic()
    inv_sub_res1 = sub_res1.inverse_mod(res)
    param = -(inv_sub_res1*sub_res0)%res
    

    return (res,param)
    
"""
Тест
"""
#K = GF(next_prime(100000000009001))
K = GF(9001)
R.<x> = PolynomialRing(K)
f = x^5+123*x^4+54*x^3+6423*x^2+625*x+467
l = 3
I = ell_torsin_ideal(f,l)



import itertools
import thermodynamic as thermo

'''
Number to return
'''
NUMBER_RETURN = 10

'''
[lenght_f, lenght_r, GC_f, GC_f, Tm_f, Tm_f, 0, 0, 0, 0, 0, 0]
'''
sc_ideal = [20, 20, 50, 50, 57, 57, 0, 0, 0, 0, 0, 0]

''' 
Weight 
'''
k = [0.5, 0.5, 1, 1, 1, 1, 0.1, 0.1, 0.2, 0.2, 0.1, 0.2]

'''
Complement
'''  
def complement(seq):
    complement = ''
    for base in seq:
        complement += complement_base(base)
    return complement


def complement_base(base):
    if (base == 'A'):
        return 'T'
    elif (base == 'T'):
        return 'A'
    elif (base == 'G'):
        return 'C'
    elif (base == 'C'):
        return 'G'
    else:
        return base
'''
Reverse
'''
def reverse(seq):
    return seq[::-1]

'''
Score of a pair of character 
'''
def s(s):
    if (s == {'A','T'}):
        return 2
    elif (s == {'C','G'}):  
        return 4
    else:
        return 0

'''
Real-valued function 
'''
def S(x, y):
    n = len(x)
    m = len(y)
    result = 0
    for k in range(-(n - 1), m):
        sum = 0
        for i in range(1, n + 1):
            j = i + k
            if (j > 0 and j < (m + 1)):
                sum += s({x[i - 1], y[j - 1]})
        if (result < sum): 
            result = sum
    return result

'''
Real-valued function 
'''
def _S(x, y):
    n = len(x)
    m = len(y)
    result = 0
    for k in range(-(n - 1), m):
        sum = 0
        for i in range(1, n + 1):
            j = i + k
            if (j >= i and j < (m + 1)):
                score = s({x[i - 1], y[j - 1]})
                if (score == 0):
                    break 
                sum += score
        if (result < sum): 
            result = sum
    return result

'''
Self annealling
'''
def sa(_p, p):
    return S(_p, p)

'''
Self-end annealling
'''
def sea(_p, p):
    return _S(_p, p)

'''
Pair annealling
'''
def pa(_p, q):
    return S(_p, q)

'''
Pair end annealling
'''
def pea(_p, q):
    return _S(_p, q)
    
'''
Weighted distance
'''
def weighted_distance(sc, sc_ideal, k):
    result = 0
    for i in range(12):
        result += k[i] * abs(sc[i] - sc_ideal[i]) 
    return result 

'''
Get primers
'''
def get_primers(seq, param, r1, r2):
    primers = set()
    forward = set()
    reverses = set()
    comp = complement(seq)
    #print('Encontrando os iniciadores F...')
    for (i, l) in r1: #candidates[0] 
        f = seq[i:i + l]
        _f = reverse(f)
        gc = round(thermo.GC(f), 2)
        tm = round(thermo.Tm1(f), 2)
        sa_ = sa(_f, f)      
        sea_ = sea(_f, f)
        if (f == 'ACTACTCAATACCGTGTACC'):
            print('GC = ' + str(gc) + ' Tm = ' + str(tm) + ' sa = ' + str(sa_) + ' sea = ' + str(sea_))
                    
        #print('Validando os iniciadores F...')
        if (param[2] <= gc and gc <= param[3]) and (param[4] <= tm and tm <= param[5]) and (sa_ <= param[6] and sea_ <= param[7]):
            forward.add((f, i, l, gc, tm, sa_, sea_))
    #print(forward)
    for (i, l) in r2: #candidates[1]
        r = comp[i:i + l]
        _r = reverse(r)
        gc = round(thermo.GC(r), 2)
        tm = round(thermo.Tm1(r), 2)
        sa_ = sa(_r, r)      
        sea_ = sea(_r, r)
        if (r == 'gagaagcaagaggagtgagt'):
            print('GC = ' + str(gc) + ' Tm = ' + str(tm) + ' sa = ' + str(sa_) + ' sea = ' + str(sea_))
        #print('Validando os iniciadores R...')
        if (param[2] <= gc and gc <= param[3]) and (param[4] <= tm and tm <= param[5]) and (sa_ <= param[6] and sea_ <= param[7]):
            reverses.add((r, i, l, gc, tm, sa_, sea_))
    #print(reverses)
    for (f, r) in itertools.product(forward, reverses):
        i_f = f[1]
        l_f = f[2]
        i_r = r[1]
        l_t = param[10]
        if (i_r - (i_f + l_f) >= l_t):
            _f = reverse(f[0])
            _r = r[0]
            pa_ = pa(_f, _r)
            pea_ = pea(_f, _r)
            #print(str(pa_) + ',' + str(pea_))
            if (pa_ <= param[8] and pea_ <= param[9]):
                primers.add((f, r, pa_, pea_))
    return primers    


'''
Key
'''    
def getKey(_list):
    return _list[0]

'''
 Main procedure
'''
def design_primer(seq, param, r1, r2):
    result = []
    # f = r = (seq, i, l, gc, tm, sa, sea)
    primers = get_primers(seq, param, r1, r2)
    for (f, r, pa, pea) in primers:
        sc = [f[2], r[2], f[3], r[3], f[4], r[4], f[5], r[5], f[6], r[6], pa, pea]
        result.append([weighted_distance(sc, sc_ideal, k), (f, r)])
    sorted_result = sorted(result, key=getKey, reverse=True)
    
    return sorted_result[0:NUMBER_RETURN]

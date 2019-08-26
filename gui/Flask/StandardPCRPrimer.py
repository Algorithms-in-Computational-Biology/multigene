import bioinformatic as bio 

'''
Complexity time: O(n^2) ??
'''
def find(l_min, l_max, begin, end):
    r = []
    for i in range(l_min, l_max + 1):
        for j in range(begin, (end + 1) - i): 
            r.append((j, i))
    return r


def design(seq, param): 
    n = len(seq)
    l_min = param[0]
    l_max = param[1]
    l_t = param[10] 

    begin = 0
    end = n - (l_t + l_min)
    forward = find(l_min, l_max, begin, end)
    begin = (l_t + l_min) 
    end = n
    reverse = find(l_min, l_max, begin, end)   
    
    return bio.design_primer(seq, param, forward, reverse)
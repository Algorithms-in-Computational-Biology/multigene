import bioinformatic as bio
from suffix import SuffixTree

''' Faz a concatenacao do conjunto de cadeias usando o separador | '''
def concat(s):
    result = ''
    m = len(s)
    for i in range(m - 1):
        result += s[i] + '|'
    result += s[m - 1]
    return result

'''
Complexity time: O(n^2) ??
'''
def find(l_min, l_max, begin, end, seq, tree):
    r = []
    for i in range(l_min, l_max + 1):
        for j in range(begin, (end + 1) - i): 
            if not tree.has_substring(seq[j:j + i]):
                r.append((j, i))
    return r   
    

def design(family, param):
    result = []
    for seq in family:
        remain = set(family)
        remain.remove(seq)
        #build GST s1|s2...|sm
        tree = SuffixTree(concat(list(remain)))
        
        n = len(seq)
        l_min = param[0]
        l_max = param[1]
        l_t = param[10] 
        
        begin = 0 
        end = n - (l_t + l_min)
        forward = find(l_min, l_max, begin, end, seq, tree)
        
        begin = (l_t + l_min) 
        end = n
        reverse = find(l_min, l_max, begin, end, seq, tree)
              
        result.append(bio.design_primer(seq, param, forward, reverse))
                
    return result

#originally shota takano 2025
import numpy as np
import re, os

### Coefficient for normalizing
std = 0.0209  #Si


def get_kpoints(option):

    ### Read PlatQlat.chk
    qlens, normconst = get_vectors()
    if qlens is None: return None

    ### Get intial k-mesh and k-space volume from option
    k_points = [float(x) * normconst for x in option]
    V = k_points[0] * k_points[1] * k_points[2]

    ### Change k-mesh according to the crystal structure (float)
    k_points = decide_k(k_points, qlens)

    ### Set the k-space volume equal to the initial value (float -> int)
    rat = np.cbrt( V / (k_points[0] * k_points[1] * k_points[2]) )
    k_points = [round(rat * x) for x in k_points]

    ### Set minimum to 3
    k_points = [max(3,k) for k in k_points]
    
    return k_points


def get_vectors():
    
    qvec, length = [], []
    with open("PlatQlat.chk","r") as f:
        txt = f.readlines()
        for i in range(0,3):
            try:
                res = re.findall(r"[-+]?\d*\.\d+|\d+", txt[i])
                qvec.append(np.array(res[3:]).astype(float))
                length.append(round(np.sqrt(sum([j*j for j in qvec[-1]])), 4))
            except:
                return None, None
                
    BZV = np.inner(qvec[0], np.cross(qvec[1], qvec[2]))
    normconst = np.cbrt( BZV / std )
    return length, normconst


def decide_k(k_points, qlens):

    ### Sort k-points and QLAT according to the size of QLAT
    sorted_indices = np.argsort(qlens)
    sorted_q = [qlens[i] for i in sorted_indices]
    sorted_k = [k_points[i] for i in sorted_indices]

    ### Cahnge k-size according to the ratio of "QLAT[i] to max(QLAT)"
    new_sorted_k = []
    for i,k in enumerate(sorted_k):
        k *= (sorted_q[i] / sorted_q[0])
        new_sorted_k.append(k)

    ### Resort k-points to original
    new_k = [0] * len(k_points)
    for i, k in enumerate(new_sorted_k):
        original_index = sorted_indices[i]
        new_k[original_index] = k

    return new_k


def get_q(k_points, kratio):
    import math
    q1 = [int(k) for k in k_points]
    q_points = [math.ceil(round(x * kratio, 3)) for x in q1]
    q_points = [max(3,q) for q in q_points]
    return q_points


### old
def decide_k0(k_points, qlens, maxIndex, medIndex, minIndex):

    ### Ratio of "max to min", "med to min"
    ratio1 = max(qlens) // min(qlens)
    ratio2 = qlens[medIndex] // min(qlens) if medIndex else None

    #[a,a,a]
    if set(maxIndex) == set(minIndex): pass  #[8,8,8]

    #[a,a,b]
    elif len(maxIndex) == 2:
       b = minIndex[0]
       if ratio1 == 1: pass  #[8,8,8]
       elif ratio1 == 2: k_points[b] = k_points[b] / 2  #[8,8,4]
       else: k_points[b] = k_points[b] / 4  #[8,8,2]

    #[a,b,b]
    elif len(minIndex) == 2:
       b1 = minIndex[0]
       b2 = minIndex[1]
       if ratio1 == 1: pass  #[8,8,8]
       elif ratio1 == 2: k_points[b1] = k_points[b2] = k_points[b1] / 2  #[8,4,4]
       else: k_points[b1] = k_points[b2] = k_points[b1] / 4  #[8,2,2]

    #[a,b,c]
    else:
       a = maxIndex[0]
       b = medIndex[0]
       c = minIndex[0]
       if ratio1 == 1: pass  #[8,8,8]
       elif ratio1 == 2:
           k_points[c] = k_points[c] / 2
           if ratio2 == 2: k_points[b] = k_points[b] / 2  #[8,4,4]
           else: pass  #[8,8,4]
       else:
           k_points[c] = k_points[c] / 4
           if ratio2 == 1: pass  #[8,8,2]
           elif ratio2 == 2: k_points[b] = k_points[b] / 2  #[8,4,2]
           else: k_points[b] = k_points[b] / 4  #[8,2,2]

    return k_points



if __name__=='__main__':
    path = "/home/stakano/result_correct/948/chk_err/"
    for mp in os.listdir(path):
        if not "mp-" in mp: continue
        os.chdir(path+mp)
        kkk = get_kpoints(option=[10,10,10])
        print(mp, kkk, get_q(kkk))


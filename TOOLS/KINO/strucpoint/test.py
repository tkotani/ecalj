import sys
from types import *

def get_var(tok):
        #input a ( b + i ) d aa -> a ( b + i )
        #input a ( b ( i ) ) f g  -> a ( b (i) )
        labels=[]
        ic=0
        icmax=ic
        for nn in tok:
                labels.append(nn)
                if nn=='(':
                        ic=ic+1
                        icmax=ic
                elif nn==')':
                        ic=ic-1
                        if ic<=0:
                                return [icmax,labels]
        # error, but continue
        return [icmax,labels]

a = ['a','(','oi','(','j',')',')','=',0]
print get_var(a)

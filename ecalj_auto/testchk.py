#!/usr/bin/env python3
import sys, os, re
from pathlib import Path

errqpu=1.1e-2
class TestChk:
    def __init__(self, previous_file, present_file, sw_print):
        "dqpu computs the differences of QPU.    requires: compareqpu"
        "dnum computs the differences of eigen values.   requires:  diffnum, linewithkey"
        self.file1 = previous_file
        self.file2 = present_file
        self.rfile1 = previous_file.read_text().split("\n")
        self.rfile2 = present_file.read_text().split("\n")
        self.min_len = min(len(self.rfile1), len(self.rfile2))
        self.sw_print = sw_print
    def dqpu(self):
        #print (' Check  QPU  of ', self.file1, 'with reference to ' ,self.file2)
        errmax = self.compareqpu()
        if errmax> errqpu:
            print( '   Error! Difference >1.1e-2 between:   ', self.file1, ' & ', self.file2, ' :  max(abs(QPU-QPU))=', errmax)
            print( '   (But small diff can be not a problem. It can be due to a minor machine-depenence.)')
            return ' QP ERR! '
        elif self.sw_print:
            print (f'   OK! Difference sum= {errmax} < {errqpu} {self.file1}')
            #print( f'       {self.file1}')
            #print( f'       {self.file2}')
        return ' QP OK! '
    def compareqpu(self):
        errmax=0.0
        sw_print = self.sw_print
        for ix in range(self.min_len):
            iline = self.rfile1[ix] if ix < len(self.rfile1) else ''
            ilin2 = self.rfile2[ix] if ix < len(self.rfile2) else ''

            if ix == 5 and sw_print:
                #print(iline)
                continue
            if ix >= 6:
                #if sw_print: print(iline[0:32],end='')
                try:
                    for iw in range(4,17):
                        w1=float(iline.split()[iw])
                        w2=float(ilin2.split()[iw])
                        diff = w1 - w2
                        #if iw <= 14 and sw_print: print( f'{diff:6.2f}', end='')
                        #elif sw_print: print( f'{diff:9.5f}', end='')
                        errmax = max(errmax, abs(diff))
                except Exception:
                    continue
        return errmax
    def dnum(self, tol, comparekeys):
        errmax = self.diffnum(tol, comparekeys)
        #print (' Check  num  of ', self.file1, ' with  reference  to ' ,self.file2)
        if errmax > tol:
            print( '   Error! MaxDiff=', errmax, '> tol=', tol, ' between   ', self.file1, ' & ', self.file2 )
            print( '   (But small diff can be not a problem. It can be due to a minor machine-depenence.)')
            return ' log ERR! '
        elif self.sw_print :
            print( f'   OK! MaxDiff={errmax}< tol={tol} {self.file1}')
            #print( f'       {self.file1}')
            #print( f'       {self.file2}')
        return ' log OK! '
    def diffnum(self, tol, comparekeys):
        errmax, ierrl = 0, 0
        sw_print = self.sw_print
        if comparekeys:
            rfile1 = self.linewithkey(comparekeys+['written','comparing'], self.rfile1)
            rfile2 = self.linewithkey(comparekeys+['written','comparing'], self.rfile2)
        for ix in range(self.min_len):
            try:
                iline = rfile1[ix].replace('D+', 'e+').replace('D-', 'e-')
                iii1 = iline.split()
                ilin2 = rfile2[ix].replace('D+', 'e+').replace('D-', 'e-')
                iii2 = ilin2.split()
            except IndexError:
                break    
            for i in range(len(iii1)):
                try:
                    w1 = float(iii1[i])
                    w2 = float(iii2[i])
                    diff = abs(w1 - w2)
                    errmax = max(errmax, diff)
                    if diff > tol:
                        if sw_print:
                            print()
                            print( ix, iline)
                            print( ix, ilin2)
                        ierrl += 1
                        if ierrl > 10:
                            print('   Error lines (less than ten shown) with > tol=', tol)
                            return errmax
                except (ValueError, IndexError):
                    continue
        return errmax
    @staticmethod
    def linewithkey(comparekeys, lines):
        matches = []
        for line in lines:
            if any(re.match(key, line) for key in comparekeys):
                matches.append(line)
        return matches            

def chk(path1, path2, sw_print):
    path1, path2 = Path(path1), Path(path2) #  print(' Compare QPU and log under ',str(path1),str(path2))
    result= ''
    d1 = path1 
    d2 = path2
    sw_print=True
    for QP in ['QPU','QPD']:
        qp1 = path1/QP   #next((f for f in d1_itr.iterdir() if re.match(r'QPU.\S+', f.name)), None)
        qp2 = path2/QP   #next((f for f in d2_itr.iterdir() if re.match(r'QPU.\S+', f.name)), None)
        if qp1.exists()!= qp2.exists():
            print('Skip! We have only one of ',qp1,qp2)
            result += ' Skip '
            continue
        elif not qp1.exists():
            continue
        testqp = TestChk(qp1, qp2, sw_print)
        result += testqp.dqpu()
        #print(QP,end='')
    try:
        log1 = path1 / next((f.name for f in d1.iterdir() if re.match(r'log\.\S+', f.name) and not f.name.startswith('log.t')), None)
        log2 = path2 / next((f.name for f in d2.iterdir() if re.match(r'log\.\S+', f.name) and not f.name.startswith('log.t')), None)
        #print(' test log '+log1.name,end='')
        testlog   = TestChk(log1, log2, sw_print)
        resultlog = testlog.dnum(tol=3e-3,comparekeys=['fp evl'])
        return result+resultlog
    except:
        result += ' Skip '
        return result

if __name__ == "__main__":
    sw_print = False
    if len(sys.argv) != 3:
        print("   usage: python testchk.py <TargetPath> <RefPath>")
        print("          Check whether <TargetPath>/mp-*/Files agree with <RefPath>/mp-*/Files.  Files are QPU, QPD and log.*")
        sys.exit(-1)
    d1root= Path(sys.argv[1])
    d2root= Path(sys.argv[2])
    directory_names = set() #directories including QPU...
    for dir_path in d1root.rglob('*'):
        if dir_path.is_dir() and re.match(r'^mp', dir_path.name): #search mp-*
            for file_path in dir_path.rglob('*'):
                if file_path.is_file() and re.search(r'^(QPU|QPD|log\.mp\S*)$', file_path.name):
                    directory_names.add(file_path.parent) #file path
    rrr=''
    if(len(directory_names)==0):
        print('ERR Nothig to compare! Check supplied directory!')
        sys.exit()
    for d1 in sorted(directory_names):
        print('---',d1, '--------------------')
        d2 = d2root / d1.relative_to(d1root)
        result = chk(d1,d2, sw_print)
        rrr+=d1.name+result+'\n'
        result = result.split()
        #if len(result)==0: print('No test performed for '+d1.name)
    #print(rrr)
    if 'ERR' in rrr:
        print('=== FAILED at some tests ===')
    elif 'Skip' in rrr:
        print('OK! ALL PASSED except Skipped!')
    else:
        print('OK! ALL PASSED!')


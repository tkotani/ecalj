class Aw_var:
        name=''
        level=0
        nameorig=''
        namesize=''
        type=''
        prevar=''
        execfuncname=''
        def __init__(self,name,type,execfuncname, level=0,nameorig='',namesize='',prevar=''):
                self.name=name
                self.type=type
                self.execfuncname=execfuncname
                self.set1(level,nameorig,namesize,prevar)


        def set1(self,level,nameorig,namesize,prevar):
                self.level=level
                self.nameorig=nameorig
                self.namesize=namesize
                self.prevar=prevar

        def __str__(self,s=','):
                n=''
                n=self.name+s+self.type+s+self.execfuncname
                n=n+s+str(self.level)+s+self.nameorig+s+self.namesize+s+self.prevar
                return n


class Awvar:
        name=''
        level=0
        nameorig=''
        namesize=''
        type=''
        prevar=''
        execfuncname=''
        def __init__(self,name,type,execfuncname,level=0,nameorig='',namesize='',prevar=''):
                self.name=name
                self.type=type
                self.execfuncname=execfuncname
                self.set1(level,nameorig,namesize,prevar)

        def set1(self,level,nameorig,namesize,prevar):
                self.level=level
                self.nameorig=nameorig
                self.namesize=namesize
                self.prevar=prevar

        def __str__(self):
                return self.name+','+str(self.level)+','+self.nameorig+','+self.namesize+','+self.type+','+self.prevar+','+self.execfuncname


a = Awvar('name','type','fname')
print a


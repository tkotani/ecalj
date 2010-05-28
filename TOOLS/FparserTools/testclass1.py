#!/usr/bin/env python


class A(object):
#     self.x='a'
     def __init__(self):
         self.x='111'
         return

class B(A):
#     self.x='b'
#     def __init__(self):
#         return
     def test3(self):
         qqq= self.test2()
         return qqq
     def test2(self):
         return 'bbbbaaaa'
     
class C(B):
     def test1(self):
         self.x='c'
         return 'c1'
#     def __init__(self):
#         self.x='333'
         return

class D(C):
     #self.x='d'
#     def __init__(self):
#         self.x='d'
     def test2(self):
         return 'd1ddd'

a=D()
print a.x
b= a.test3()
print b


c=C()
print c.test3()

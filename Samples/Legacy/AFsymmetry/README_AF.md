We can perform LDA/QSGW calculation with the AF symmetry.
Only up spin are calculated. Then charge density and eigenvalues are
symmetrized for the antiferro symmetry.


To set AF symmetry,
We have to set a line in ctrl file as
```
SYMGRPAF i:(1,1,1)  #Antiferro symmetry operation
```
Here SYMGRPAF is the antiferro magnetic symmetry that 
`g =  i:(1,1,1) * i $\sigma$_y$ should be the generator of magnetic space group,
whereas 
```
SITE    ATOM=Niup POS=  .0   .0   .0  AF=1    #Antiferro pair
        ATOM=Nidn POS= 1.0  1.0  1.0  AF=-1
```
should contains the `AF=` index to specify antiferro pairs.





(Caution. I do not think SYMGRPAF works well when SOC is fully included, so=1)



This is an example of ctrl.nise. NiSe. 
```
  SYMGRPAF  i:(0,0,1/2)
  ATOM=Niup POS=     0.00       0.000   0.000 AF=1
  ATOM=Nidn POS=     0.00       0.000   0.500 AF=-1
```
, where inversion with (0,0,1/2) translation gives AF symmetry.








import os
for i in os.listdir('.'):
    if not ('list.' in i): continue
    print(i)
    ccc='sed -e s/mmmmm/'+i+'/g jobtemplate >jobx.'+i
    print(ccc)
    os.system(ccc)
#qsub

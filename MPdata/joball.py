import os
for i in os.listdir('.'):
    if not ('list.' in i): continue
    if ('jobx' in i): continue
    if ('out' in i): continue
    if ('log' in i): continue
    if ('atom' in i): continue
    print(i)
    ccc='sed -e s/mmmmm/'+i+'/g jobtemplate >jobx.'+i
    print(ccc)
    os.system(ccc)
    os.system('qsub jobx.'+i)
#qsub

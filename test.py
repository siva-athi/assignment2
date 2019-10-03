import re

from numpy.ma import count

seq='aaattatagggatatata'

motif='ata'

Q=re.compile(motif)

a='('+'\''+motif+'\','+str(count([item.start(0) for item in Q.finditer(seq)]))+')'
print(a)

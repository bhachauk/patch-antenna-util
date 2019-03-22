import os
import sys

scriptpath = "patch.py"

sys.path.append(os.path.abspath(scriptpath))

# Do the import
from patch import DesignPatch, inputImpedance, insetFeedPosition, getDirectivity, PatchEHPlanePlot, SurfacePlot, getGs


freq = 2.4e9
Er = 4.4                                                          
h = 1.6* 10 ** -3
v = 3 * 10 ** 8

W, L = DesignPatch(Er, h, freq)


Rin = inputImpedance(freq, W, L, h, Er)  

print 'Inset Feed Position : ', insetFeedPosition(Rin, L)

G1, G12 = getGs(freq, W, L)

print 'G1 : ', G1

print 'G12 : ', G12

I1=1.863

I2=3.59801

d1, d2 = getDirectivity(G1, G12, W, freq, I1 , I2)                            

print 'Directivity : ', d1, ' dB'
print 'Directivity : ', d2, ' dB'



fields = PatchEHPlanePlot(freq, W, L, h, Er)


SurfacePlot(fields, freq, W, L, h, Er)
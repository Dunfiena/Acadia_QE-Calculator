#!/usr/bin/env python
from bands import *

datafile=
fermi =
symmetryfile=
bool_shift_efermi= True
fig, ax = plt.subplots()

#bndplot(datafile,fermi,symmetryfile,ax)
bndplot(datafile,fermi,symmetryfile,ax,shift_fermi=0,\
color='black',linestyle='solid',name_k_points=['Γ','X','W','K','Γ','L','U','W','L','K|U','X'],legend='Si, PBE')


fig.savefig


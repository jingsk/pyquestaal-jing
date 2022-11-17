from ase.io import read
import numpy as np

#si = read('POSCAR')
#kpath = si.cell.bandpath('GXWG', npoints=100)
def write_syml(kpath, line=False, name='temp'):
    if not line:
        head='#written from path='+kpath.path+'with nkpt='+str(int(kpath.kpts.shape[0]))
        np.savetxt("syml."+name,np.array(kpath.kpts), fmt='%.4f', header=head)
    else:
        f = open("syml."+name, 'w')
        four_dec={'float_kind':lambda x: "%.4f" % x}
        print('#written from path=',kpath.path,'with nkpt=',str(int(kpath.kpts.shape[0])),file=f)
        inds=sorted(kpath._find_special_point_indices(eps=1e-5).keys())
        for i in range(len(inds)-1):
            nkpt=inds[i+1]-inds[i]+1
            kpt1=np.array2string(np.array(kpath.kpts[inds[i]]), formatter=four_dec)
            kpt2=np.array2string(np.array(kpath.kpts[inds[i+1]]), formatter=four_dec)
            print(nkpt,"   ",str(kpt1).strip("[]"),"   ",str(kpt2).strip("[]"),"   ",kpath.path[i]," to ",kpath.path[i+1],file=f)
        f.close()
    return
#write_syml(kpath,name='si')

#!/usr/bin/env python

import numpy as np

if __name__ == '__main__':
    import sys
    import datetime

    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print("")
    print(">>>>> Code by Young Jae Choi @ POSTECH <<<<<".center(120))
    print(("code started time: "+time).center(120))
    print("")
    print("==================================================================================================".center(120))
    print("")
    print("Useage  ==> ./ase-cell-magnifier.py 'structure file' 'Mxx' 'Mxy' ..(omitted).. 'Mzy' 'Mzz'".center(120))
    print("Example ==> ./ase-cell-magnifier.py 'structure file' 1 0 0   0 0.9 0.1   0 0.1 0.9".center(120))
    print("             ( Mij is transformation matrix; A_new = M @ A_old )                                ".center(120))
    print("    Or  ==> ./ase-cell-magnifier.py 'structure file' 'Mxx' 'Myy' 'Mzz'".center(120))
    print("Example ==> ./ase-cell-magnifier.py 'structure file'  1.2   1.1   1".center(120))
    print(" >>>>     Note) 'SCALED_position' remains unchanged     <<<<".center(120))
    print("")
    print("==================================================================================================".center(120))
    
    #### Define transformation matrix
    if len(sys.argv) is 11:
        M = np.reshape([float(sys.argv[i+2]) for i in range(9)], (3,3))
    elif len(sys.argv) is 5:
        M = np.zeros((3,3))
        for i in range(3):
            M[i,i] = sys.argv[i+2]
    else:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(120))
        sys.exit(1)
    print("")

    #### Load structure
    a_file = sys.argv[1]
    from ase.io import read, write
    alist = read(a_file, index=':')
    if not isinstance(alist, list):
        alist = [alist]
    
    for i in range(len(alist)):
        atoms = alist[i]
        A     = atoms.get_cell()
        A_new = np.matmul(M, A)
        atoms.set_cell(A_new, scale_atoms=True)

    if a_file[-5:] == '.traj':
        o_file = 'magnified_'+a_file
    else:
        o_file = 'magnified_'+a_file+'.traj'
    write(o_file, alist)


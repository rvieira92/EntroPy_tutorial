import numpy as _np
def ReadTsweep(fileName):
    """
    #Temperature   <M>   <M^2>   <M^4>   U_{Binder}   \chi   C_v(tot)   <E> ...
               

    Args:
        fileName (_type_): _description_

    Returns:
        _type_: _description_
    """
    data = _np.genfromtxt(fileName,usecols=[0,1,4,5,7])
    # 0, T : 1, <M> (mu_b) : 4, U_Binder : 5, chi : 7, E (mRy/atom)
    data = data[np.argsort(data[:,0])]
    c = np.gradient(data[:,-1]/(1E3*kb_Ry),data[:,0],axis=0)
    data=np.c_[data,c]

    return data
       
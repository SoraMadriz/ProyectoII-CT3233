import numpy as np #Dependencia de la estructura

def functionYbus(super_matrix,num_bar,num_nexos):
    a_pr = np.zeros((num_nexos,num_bar))        #Mat. incidencia

#CONSTRUCCION DE LA MATRIZ DE ADMITANCIAS  
    for nexus,arista in enumerate(super_matrix[:,:-1].tolist()):
        #vertices de entrada y salida
        vertex_out = int(arista[0] -1)
        vertex_in = int((arista[1]-1) if arista[1]!=0 else arista[1])
        #actualización de la matriz de incidencia
        a_pr[nexus,vertex_in] = ((-1) if vertex_in!=0 else 0)
        a_pr[nexus,vertex_out] = 1

#CONSTRUCCION DE LA MATRIZ IMPEDANCIA DE RAMA --> ADMITANCIA DE RAMA
    y_rama = np.diagflat(super_matrix[:,-1]).astype(complex)

#CALCULO DE YBUS
    a_pr_trans = np.transpose(a_pr)
    y_bus = (a_pr_trans)@(y_rama)@(a_pr)

    return y_bus
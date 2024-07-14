import numpy as np
import pdb

def linealizared(p_gen,p_load,z_lines,z_trx,tap,bar_type):
#VALORES ESPECIFICOS
    p_i = p_gen - p_load

#REACTANCIAS DE LINEAS Y TRX
    if len(z_trx) != 0:
        x_trx = z_trx[:,1:]
        x_trx[:,-1] = [x/t for x,t in zip(x_trx[:,-1],tap)]
        x_total = np.concatenate((z_lines,x_trx))
    else:
         x_total = z_lines
 
#CONSTRUCCION DE LA MATRIZ DE REACTANCIAS
    # pdb.set_trace()
    react_matrix = np.zeros((len(p_i),len(p_i)))
    #Fuera de la diagonal
    for k in range(len(x_total)):
        i = x_total[k,0] -1; j = x_total[k,1]-1
        react_matrix[i,j] = -1/x_total[k,2]; react_matrix[j,i] = react_matrix[i,j]
    #Valores de la diagonal
    react_matrix = react_matrix - np.diag(np.sum(react_matrix,axis=1))

#ELIMINAR BARRAS SLACKS
    #Eliminando barras SLACK
    for i in range(len(z_lines)):
        if bar_type[i] == "SL":
            react_matrix = np.delete(react_matrix,i,0)
            react_matrix = np.delete(react_matrix,i,1)
            p_i = np.delete(p_i,i)

#OBTENER ANGULO DE LAS BARRA 
    bar_angle = (np.linalg.inv(react_matrix)) @ (np.transpose(p_i))
    localetion = np.where(bar_type =="SL")[0]
    bar_angle = np.insert(bar_angle,localetion,0)
    bar_angle = bar_angle*180/np.pi

    return bar_angle
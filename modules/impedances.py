from modules.readdata import *

#Paralelo de elementos
def parallels(y):
    #Ordenar matriz de mayor a menor
    y_ordenada = np.concatenate(y[np.argsort(y[:, 1],axis=0,kind='heapsort')],axis=0)
    i_values = []
    #Algoritmo para calcular el paralelo de admitancias
    for i in range(len(y_ordenada)):
        if (i!=0) and (y_ordenada[i,1] == y_ordenada[i-1,1]) and (y_ordenada[i,2] == y_ordenada[i-1,2]):
            new_row = y_ordenada[i,3] + y_ordenada[i-1,3]
            i_values.append(i-1)
            y_ordenada[i,3] = new_row
    y_ordenada = np.delete(y_ordenada,i_values,axis=0)

    return y_ordenada


#MODELADOS DE LAS LINEAS
z_lines = df_lines[:,(1,2,4)].copy()
y_lines = df_lines[:,:4].copy()
b_sht = df_lines[:,(1,2,-1)].copy()
y_lines[:,3] =1/(df_lines[:,3] + df_lines[:,4]*1j)
sht_lines_i = np.matrix([[row[0,0],0,row[0,-1]*1j/2] for row in b_sht],dtype=object)
sht_lines_j = np.matrix([[row[0,1],0,row[0,-1]*1j/2] for row in b_sht],dtype=object)
sht_lines = np.concatenate((sht_lines_i,sht_lines_j),axis=0)

#MODELADOS DE LAS TRX
z_trx = df_trx[:,:-2].copy()
position=np.array(df_trx[:,-1].copy()).flatten()
#Verificar posicion del tap
for k in range(len(z_trx)):
    if position[k] != z_trx[k,1]:
        z_trx[k,1], z_trx[k,2] = z_trx[k,2], z_trx[k,1]
y_trx= z_trx.copy()
tap_trx = df_trx[:,-2].copy()
trx_earth_i = np.matrix(np.zeros((len(z_trx),3),dtype=object))
trx_earth_i[:,0] = z_trx[:,1].copy()
trx_earth_j = np.matrix(np.zeros((len(z_trx),3),dtype=object))
trx_earth_j[:,0] = z_trx[:,2].copy()
for i in range(len(z_trx)):
    if tap_trx[i] == 1:
        y_trx[i,3] = (tap_trx[i,0])*1/(df_trx[i,3]*1j)
    else:
        y_trx[i,3] = (tap_trx[i,0])*1/(df_trx[i,3]*1j)

        aux1 = (1-tap_trx[i,0])*1/(df_trx[i,3]*1j)
        trx_earth_i[i,2] = aux1
        aux2 = ((tap_trx[i,0])**2-tap_trx[i,0])*1/(df_trx[i,3]*1j)
        trx_earth_j[i,2] = aux2
trx_earth = np.concatenate((trx_earth_i,trx_earth_j),axis=0)



#MODELADO DE LOS ELEMENTOS SHUNT
z_sht = df_sht[:,:3].copy()
y_sht = df_sht[:,:3].copy()
y_sht[:,2] = 1/(df_sht[:,2] + df_sht[:,3]*1j)
y_sht = np.insert(y_sht,2,0,axis=1)
if (len(y_sht) != 0) and (len(y_sht) > 1):
    y_sht = parallels(y_sht)

#VALORES CONECTADAS A TIERRA
y_earth = [arr for arr in [y_sht, trx_earth, sht_lines] if arr.size>0]
y_earth = np.concatenate(y_earth)
y_earth = np.delete(y_earth,np.where(y_earth[:,2] == 0)[0],0)

if (len(y_earth) != 0) and (len(y_earth) > 1):
    y_earth = np.insert(y_earth,0,None,axis=1)
    y_earth = parallels(y_earth)

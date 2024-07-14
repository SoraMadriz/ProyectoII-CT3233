import numpy as np
import pdb

def powerLuxury(voltages ,y_bus,i_lines,j_lines,sht_lines,y_trx,sht_trx_i, sht_trx_j):
#ARRAYS DE SALIDA
    index = list(zip(i_lines,j_lines))
    flux_ij = np.array([]); flux_ji = np.array([])
    

#FLUJO DE POTENCIA POR NEXOS - LINEAS
    for position,nexus in enumerate(index):
    #Indices posicionales
        i = nexus[0] - 1
        j = nexus[1] - 1
    #Flujo por el nexo ij
        #Valor de las capacitancias de las líneas
        y_sht_i = sht_lines[position,2]/2
        #Conexiones entre barras ij
        term1 = abs(voltages[i])**2*y_sht_i
        term2 = voltages[i]*np.conjugate((voltages[i]-voltages[j])*(-1*y_bus[i,j]))
        #Flujo ij
        flux_ij = np.append(flux_ij,term1 + term2)
    #Flujo por el nexo ji
        #Valor de las capacitancias de las líneas
        y_sht_j = sht_lines[position,2]/2
        #Conexiones entre barras ij
        term1 = abs(voltages[j])**2*y_sht_j
        term2 = voltages[j]*np.conjugate((voltages[j]-voltages[i])*(-1*y_bus[i,j]))
        #Flujo ij
        flux_ji = np.append(flux_ji,term1 + term2)

#FLUJO DE POTENCIA POR NEXOS - TRX
    if len(y_trx) != 0:
        for position,nexus in enumerate(y_trx[:,1:]): 
            nexus=np.array(nexus).flatten()
            i = nexus[0] - 1
            j = nexus[1] - 1
    #Valor de las capacitancias de las líneas
            y_i = sht_trx_i[position,2]
            y_j = sht_trx_j[position,2]
    #Flujo por el nexo ij
        #Conexiones entre barras ij
            term1 = abs(voltages[i])**2*y_i
            term2 = voltages[i]*np.conjugate((voltages[i]-voltages[j])*y_trx[position,3])
        #Flujo ij
            flux_ij = np.append(flux_ij,term1 + term2)
    #Flujo por el nexo ji
        #Conexiones entre barras ij
            term1 = abs(voltages[i])**2*y_j
            term2 = voltages[j]*np.conjugate((voltages[j]-voltages[i])*y_trx[position,3])
        #Flujo ij
            flux_ji = np.append(flux_ji,term1 + term2)

#EXTRACCION DE POTENCIA ACTIVA Y REACTIVA
    #Potencia Activa
    p_ij = np.real(flux_ij)
    p_ji = np.real(flux_ji)
    #Potencia Reactiva
    q_ij = np.imag(flux_ij)
    q_ji = np.imag(flux_ji)

#PERDIDA POR LOS NEXOS
    #Potencia activa
    p_losses = p_ij + p_ji
    q_losses = q_ij + q_ji
    
    return p_ij,p_ji,q_ij,q_ji, p_losses, q_losses
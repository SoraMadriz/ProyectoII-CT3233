import numpy as np

def gaussSeidel(tolerance, max_iter,y_bus,p_gen,p_load,q_gen,q_load,v_bar,bar_type):
#INICIALIZACION DE VARIABLES
    v_salida = np.zeros((len(y_bus)),dtype=object)
    error = 0

#VALORES ESPECIFICADOS
    p_i = p_gen - p_load
    q_i = q_gen - q_load

#METODO DE GAUSS-SEIDEL
    for k in range(max_iter):                                    
        for i in range(len(y_bus)):                                            
            if bar_type[i] == 'SL':                                 
                v_salida[i] = v_bar[i]
                
            elif bar_type[i] == 'PQ': 
                term1=0; sigma1=0; sigma2=0; sigma_q = 0                            
                #Primera sumatoria
                for j in range(0,len(v_salida)-1):
                    if i != j:
                        sigma1 += y_bus[i,j]*v_salida[j]
                    else:
                        break
                #Segunda sumatoria
                for j in range(i,len(v_bar)):
                    if i != j:
                        sigma2 += y_bus[i,j]*v_bar[j]
                term1 = np.conj((p_i[i] + 1j*q_i[i])/v_bar[i])
                v_salida[i] = 1/y_bus[i,i]*(term1 - (sigma1 + sigma2))

            elif bar_type[i] == 'PV':
                V_DATO = abs(v_bar[i]) 
                term1=0; sigma1=0; sigma2=0; sigma_q = 0
                #Correccion Q_pv
                #Primera sumatoria
                for j in range(0,len(v_salida)):
                    if i != j:
                        sigma1 += y_bus[i,j]*v_salida[j]
                        sigma_q += y_bus[i,j]*v_salida[j]
                    else:
                        break
                #Segunda sumatoria
                for j in range(i,len(v_bar)):
                    if i != j:
                        sigma2 += y_bus[i,j]*v_bar[j]
                        sigma_q += y_bus[i,j]*v_bar[j]
                    else:
                        sigma_q += y_bus[i,j]*v_bar[j]
                #Expresion de q para barras PV
                q_i[i] = -1*(np.imag(np.conj(v_bar[i])*(sigma_q)))
                #Calculo del voltaje
                term1 = np.conj((p_i[i] + 1j*q_i[i])/v_bar[i])
                v_salida[i] = (term1 - (sigma1 + sigma2))*1/y_bus[i,i]
                v_salida[i] = (v_salida[i]/abs(v_salida[i]))*V_DATO

        error = max(abs(v_salida-v_bar))
        if error < tolerance:
            break
        else:
            v_bar = v_salida.copy() 
    return v_salida, p_i, q_i,[float(k)],[error]
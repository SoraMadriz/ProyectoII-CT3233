import numpy as np
import pdb

def deliveredpower(voltages,y_bus,type_bar):
#POTENCIA GENERADA
    S_calc = np.array([],dtype=object)
    for i in range(len(y_bus)):
        # pdb.set_trace()
        if type_bar[i] == "SL":
    #Termino que se encuentra fuera de la sumatoria:
            term1 = (abs(voltages[i])**2)*np.conjugate(y_bus[i,i])
    #Sumatoria
            term2 = 0
            for j in range(len(y_bus)):
                if j != i:
                    term2 += np.conjugate(voltages[j]*y_bus[i,j])
            S_calc = np.append(S_calc, term1 + voltages[i]*term2)

        elif type_bar[i] == "PV":
    #Termino que se encuentra fuera de la sumatoria:
            term1 = (abs(voltages[i])**2)*np.conjugate(y_bus[i,i])
    #Sumatoria
            term2 = 0
            for j in range(len(y_bus)):
                if j != i:
                    term2 += np.conjugate(voltages[j]*y_bus[i,j])
            S_calc = np.append(S_calc, np.imag(term1 + voltages[i]*term2)*1j)
        elif type_bar[i] == "PQ":
            S_calc = np.append(S_calc,0+0j)
    return S_calc
    

import numpy as np
import sympy as sp

def quickDecoupling(tolerance, max_iter,y_bus,p_gen,p_load,q_gen,q_load,v_bar,bar_type):
#VALORES DE SALIDA
    error = 0
    equation_p = np.zeros((len(y_bus)),dtype=object)
    equation_q = np.zeros((len(y_bus)),dtype=object)
    values = {}
    index_delete_sl = []
    index_delete_pv = []

#VALORES ESPECIFICADOS
    p_i = p_gen - p_load; p_return = p_i
    q_i = q_gen - q_load; q_return = q_i

#SEPARANDO LA MATRIZ DE ADMITANCIA Y VOLTAJES EN SUS MODULOS Y ANGULOS
    #Matriz de admitancia
    ybus_mod = abs(y_bus)
    ybus_ang = np.angle(y_bus)
    ybus_imag = np.imag(y_bus)
    #Voltajes
    vbar_mod = sp.Matrix(abs(v_bar))
    vbar_ang = sp.Matrix(np.arctan2(np.imag(v_bar),np.real(v_bar)))

#CREACION DE LAS VARIABLES A CALCULAR
    #VOLTAJES
    vol_var = [sp.symbols(f'V{i+1}') for i in range(len(v_bar))]
    #Angulos
    ang_var = [sp.symbols(f'alpha{i+1}') for i in range(len(v_bar))]

#ECUACIONES DEL CALCULO ITERATIVO
    #Ecuaciones de potencia
    for i in range(len(bar_type)):
        if bar_type[i] == "SL":
            index_delete_sl.append(i)
        elif bar_type[i] == "PQ":
            for j in range(len(bar_type)):
                equation_p[i] += vol_var[i]*ybus_mod[i,j]*vol_var[j]*sp.cos(ybus_ang[i,j] - ang_var[i] + ang_var[j])
                equation_q[i] -= vol_var[i]*ybus_mod[i,j]*vol_var[j]*sp.sin(ybus_ang[i,j] - ang_var[i] + ang_var[j])
        elif bar_type[i] == "PV":
            for j in range(len(bar_type)):
                equation_p[i] += vol_var[i]*ybus_mod[i,j]*vol_var[j]*sp.cos(ybus_ang[i,j] - ang_var[i] + ang_var[j])
            index_delete_pv.append(i)

    #Valores a no eliminar
    no_delete_sl = set([i for i, _ in enumerate(vbar_mod) if i not in index_delete_sl])
    no_delete_pv = set([i for i, _ in enumerate(vbar_mod) if i not in index_delete_pv])

    #Añadiendo la potencia especifica a cada ecuacion de potencia
    equation_p = sp.Matrix(p_i) - sp.Matrix(equation_p)
    equation_q= sp.Matrix(q_i) - sp.Matrix(equation_q)

    #Limpieza de variables:
    equation_q = np.delete(equation_q,index_delete_sl + index_delete_pv)
    v_mod = sp.Matrix(np.delete(vbar_mod,index_delete_pv + index_delete_sl))
    v_ang = sp.Matrix(np.delete(vbar_ang,index_delete_sl))
    p_i = np.delete(p_i,index_delete_sl)
    q_i = np.delete(q_i,index_delete_sl + index_delete_pv)
    vol_var = np.delete(vol_var,index_delete_pv + index_delete_sl)
    ang_var = np.delete(ang_var, index_delete_sl)

    #Construccion de las matrices de susceptancias
    b_angle = sp.Matrix(np.delete(ybus_imag,list(index_delete_sl),axis=1))
    b_angle = sp.Matrix(np.delete(b_angle,list(index_delete_sl),axis=0))
    b_angle = b_angle.inv()

    b_mod = sp.Matrix(np.delete(ybus_imag, list(index_delete_pv+index_delete_sl),axis=1))
    b_mod = sp.Matrix(np.delete(b_mod, list(index_delete_pv+index_delete_sl),axis=0))
    b_mod = b_mod.inv()

    #Vectores de ecuaciones delta
    deltaQV = sp.Matrix([equation_q[i]/abs(v_mod[i]) for i in range(len(v_mod))])
    deltaPV = sp.Matrix([equation_p[i]/abs(v_bar[i]) for i in range(len(v_bar)) if i not in index_delete_sl])

#ALGORITMO DE DESACOPLADO RAPIDO
    for w in range(max_iter):
        for i in range(len(vbar_mod)):
            values[f'V{i+1}'] = vbar_mod[i]
            values[f'alpha{i+1}'] = vbar_ang[i]

    #Sustituimos valores
        f_eval_angle = deltaPV.subs(values)
        f_eval_mod = deltaQV.subs(values)


    #Expresion
        values_next_angle = v_ang - b_angle*f_eval_angle
        values_next_mod = v_mod - b_mod*f_eval_mod

    #Comprobación de condiciones
        error_mod = max(abs(values_next_mod - v_mod))
        error_ang = max(abs(values_next_angle - v_ang))
        error = max(error_mod,error_ang)

        if error < tolerance:
            break
        else:
            k=0
            for i in list(no_delete_sl.intersection(no_delete_pv)):
                vbar_mod[i] = values_next_mod[k]
                k += 1
    
            k=0 #iterador auxiliar
            for i in list(no_delete_sl):
                vbar_ang[i] = values_next_angle[k]
                k += 1

            v_ang = values_next_angle.copy()
            v_mod = values_next_mod.copy()
    return vbar_mod, vbar_ang, p_return, q_return, [w], [error]
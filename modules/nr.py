import numpy as np
import sympy as sp
import pdb  
def newtonRaphson(tolerance, max_iter,y_bus,p_gen,p_load,q_gen,q_load,v_bar,bar_type):

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
    #Voltajes
    vbar_mod = sp.Matrix(abs(v_bar))
    vbar_ang = sp.Matrix(np.arctan2(np.imag(v_bar),np.real(v_bar)))
    values_bf = vbar_ang.col_join(vbar_mod)

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

    #Limpieza de variables:
    equation_p = np.delete(equation_p,index_delete_sl)
    equation_q = np.delete(equation_q,index_delete_sl + index_delete_pv)
    v_mod = sp.Matrix(np.delete(vbar_mod,index_delete_pv + index_delete_sl))
    v_ang = sp.Matrix(np.delete(vbar_ang,index_delete_sl))
    p_i = np.delete(p_i,index_delete_sl)
    q_i = np.delete(q_i,index_delete_sl + index_delete_pv)
    vol_var = np.delete(vol_var,index_delete_pv + index_delete_sl)
    ang_var = np.delete(ang_var, index_delete_sl)

    #Añadiendo la potencia especifica a cada ecuacion de potencia
    equation_p = sp.Matrix(p_i) - sp.Matrix(equation_p)
    equation_q= sp.Matrix(q_i) - sp.Matrix(equation_q)


    #Jacobianos
    #jacobianos de potencia activa
    JdeltaPalpha = equation_p.jacobian([values for values in ang_var])
    JdeltaPvolt = equation_p.jacobian([values for values in vol_var]) 
    jacobian_p = JdeltaPalpha.row_join(JdeltaPvolt)
    #Jacobiano de potencia reactiva
    JdeltaQalpha = equation_q.jacobian([values for values in ang_var])
    JdeltaQvolt = equation_q.jacobian([values for values in vol_var])
    jacobian_q = JdeltaQalpha.row_join(JdeltaQvolt)

    #Uniendo valores
    values_bf = v_ang.col_join(v_mod)
    Jacobian = jacobian_p.col_join(jacobian_q)
    f_powers = equation_p.col_join(equation_q)
    values_next = sp.zeros(values_bf.rows,values_bf.cols)


#ALGORITMO DE NR
    for k in range(max_iter):
        for i in range(len(vbar_mod)):
            values[f'V{i+1}'] = vbar_mod[i]
            values[f'alpha{i+1}'] = vbar_ang[i]

    #Sustituimos valores
        Jacobian = Jacobian.subs(values)
        Jacobian_eval = Jacobian.inv()
        f_eval = f_powers.subs(values)

    #Expresion de NR
        values_next = values_bf - Jacobian_eval*f_eval

        error = max(abs(values_next - values_bf))

    #Actualizacion de variables
        #modulos de voltajes
        k=1 #iterador auxiliar
        for i in list(no_delete_sl.intersection(no_delete_pv)):
            vbar_mod[i] = values_next[(len(v_ang)-1)+k]
            k += 1
    
        k=0 #iterador auxiliar
        for i in list(no_delete_sl):
            vbar_ang[i] = values_next[k]
            k += 1
        values_bf = values_next.copy()
        
        if error < tolerance:
            break

    return vbar_mod, vbar_ang, p_return, q_return, [k], [error]
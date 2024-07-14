from modules.impedances import *
from modules.ybus import functionYbus as fbus
from modules.gs import gaussSeidel as fgs
from modules.sflow import powerLuxury as fnexus
from modules.nr import newtonRaphson as fnr
from modules.fd import quickDecoupling as fqd
from modules.dc import linealizared as fline
from modules.sgbal import deliveredpower as fpower
from modules.writedata import writenResults as write

def main():
#-----------------------------VARIABLES DEL ARCHIVO READDDATA E IMPEDANCES---------------------------------
    #Configuraciones
    converge_radius = df_config[5,1]                #Tolerancia del programa    
    iterations_max = df_config[6,1]                 #Máxima cantidad de iteraciones permitidas
    option = np.array(df_config[:,-1]).flatten()    #Selección de los métodos a ejecutar
    option = np.delete(option,[4,5,6,7,8,9,10])     #Eliminar valores que no nos sirve para la selección
    SALIDA = df_config[-1,-1]
    print(option)                 

    #Bus
    bar_type = np.squeeze(np.asarray(df_bus[:,2]))              #Clasf.Barras del sistema
    p_load = np.squeeze(np.asarray(df_bus[:,6]))                #Potencia activa demandada
    q_load = np.squeeze(np.asarray(df_bus[:,7]))                #Potencia reactiva demandada
    p_gen = np.squeeze(np.asarray(df_bus[:,4]))                 #Potencia activa generada
    q_gen = np.squeeze(np.asarray(df_bus[:,5]))                 #Potencia reactiva generada
    v_bar = np.squeeze(np.asarray(df_bus[:,3],dtype=complex))   #Voltajes en la barra

    #Lines
    index_i_lines = np.array(df_lines[:,1]).flatten()           #Conex.Líneas en las barras i
    index_j_lines = np.array(df_lines[:,2]).flatten()           #Conex.Líneas en las barras j
    y_sht = df_lines[:,[1,2,5]]; y_sht[:,-1] = y_sht[:,-1]*1j   #Capacitancias de ls líneas

    #Trx            
    index_i_trx = np.array(df_trx[:,1]).flatten()               #Conex.TRX en las barras i
    index_j_trx = np.array(df_trx[:,2]).flatten()               #Conex.TRX en las barras j

    #Impedances
    if (len(y_trx) > 0) and (len(y_earth) > 0):
        super_matrix_yb = np.concatenate((y_lines[:,1:],y_trx[:,1:],y_earth[:,1:]))
    elif (len(y_trx) == 0) and (len(y_earth) > 0):
        super_matrix_yb = np.concatenate((y_lines[:,1:],y_earth[:,1:]))
    elif (len(z_trx) > 0) and (len(y_earth) == 0):
        super_matrix_yb = np.concatenate((y_lines[:,1:],y_trx[:,1:]))
    else:
        super_matrix_yb = np.concatenate(y_lines[:,1:])

    #Dimensionamiento de Ybus
    if len(y_trx) != 0:
        num_bar_yb = max(np.max(y_lines[:,1]),np.max(y_lines[:,2]),np.max(y_trx[:,1]),np.max(y_trx[:,2]))
    elif len(y_trx) == 0:
        num_bar_yb = max(np.max(y_lines[:,1]),np.max(y_lines[:,2]))

    num_nexos_yb = len(y_lines) + len(y_earth) + len(y_trx)

#-------------------------------------------METODOS DE ESTUDIO---------------------------------------------
    #METODO DE INCIDENCIA NODAL (YBUS)
    y_bus = fbus(super_matrix_yb,num_bar_yb,num_nexos_yb)

    #FLUJO DE POTENCIA POR EL METODO DE GAUSS-SEIDEL
    if option[0] == "Y":
        flow_gs,p_gs,q_gs,iteration,error = fgs(
            converge_radius,
            iterations_max,
            y_bus,
            p_gen,
            p_load,
            q_gen,
            q_load,
            v_bar,
            bar_type
        )
    #Valores de escritura de Gauss-Seidel
        module_vgs = [float(np.abs(i)) for i in flow_gs] 
        angle_vgs = [float(np.angle(i)*180/np.pi) for i in flow_gs]
    #Potencia entregada
        deliveried_power = fpower(flow_gs,y_bus,bar_type)
        p_gen_gs = [float(np.real(i)) for i in deliveried_power]
        q_gen_gs = [float(np.imag(i)) for i in deliveried_power]
    #Flujo de potencia por nexos
        p_ij,p_ji,q_ij,q_ji, p_losses, q_losses = fnexus(
            flow_gs, y_bus, index_i_lines, index_j_lines,
            y_sht, y_trx, trx_earth_i, trx_earth_j
            )
    #PRUEBA ----> ELIMINAR <-----
        # print("Modulo",[round(float(np.abs(i)),5) for i in flow_gs])
        # print("Ángulo",[round(float(np.angle(i)*180/np.pi),2) for i in flow_gs])
        

    #FLUJO DE POTENCIA POR EL METODO DE NEWTON-RAPHSON
    if option[1] == "Y":
        flow_nr_mod,flow_nr_ang, p_nr, q_nr, iteration_nr, error_nr = fnr(
            converge_radius,
            iterations_max,
            y_bus,
            p_gen,
            p_load,
            q_gen,
            q_load,
            v_bar,
            bar_type
        )
    #Valores de escritura de Newton-Raphson
        flow_nr_mod = [float(i) for i in np.array(flow_nr_mod).flatten()]
        flow_nr_ang = [float(i) for i in np.array(flow_nr_ang).flatten()]
        error_nr = [float(i) for i in error_nr]
        flow_nr = [mod*(np.cos(ang)+1j*np.sin(ang)) for mod,ang in zip(flow_nr_mod,flow_nr_ang)]
    #Potencia entregada
        deliveried_power = fpower(flow_nr,y_bus,bar_type)
        p_gen_nr = [np.real(i) for i in deliveried_power]
        q_gen_nr = [np.imag(i) for i in deliveried_power]
    #Flujo de potencia por nexos
        p_ij_nr,p_ji_nr,q_ij_nr,q_ji_nr, p_losses_nr, q_losses_nr = fnexus(
            flow_nr, y_bus, index_i_lines, index_j_lines,
            y_sht, y_trx, trx_earth_i, trx_earth_j
            )

    #FLUJO DE POTENCIA POR EL METODO DE DESACOPLADO RAPIDO
    if option[2] == "Y":
        flow_fd_mod, flow_fd_ang,p_fd,q_fd,iteration_fd,error_fd = fqd(
            converge_radius,
            iterations_max,
            y_bus,
            p_gen,
            p_load,
            q_gen,
            q_load,
            v_bar,
            bar_type
        )
        flow_fd_mod = [float(i) for i in np.array(flow_fd_mod).flatten()]
        flow_fd_ang = [float(i) for i in np.array(flow_fd_ang).flatten()]
        flow_fd = [mod*(np.cos(ang)+1j*np.sin(ang)) for mod,ang in zip(flow_fd_mod,flow_fd_ang)]
        error_fd = [float(i) for i in error_fd]
    #Potencia entregada
        deliveried_power = fpower(flow_fd,y_bus,bar_type)
        p_gen_fd = [np.real(i) for i in deliveried_power]
        q_gen_fd = [np.imag(i) for i in deliveried_power]
    #Flujo de potencia por nexos
        p_ij_fd,p_ji_fd,q_ij_fd,q_ji_fd, p_losses_fd, q_losses_fd = fnexus(
            flow_fd, y_bus, index_i_lines, index_j_lines,
            y_sht, y_trx, trx_earth_i, trx_earth_j
            )
        
    #METODO DE LINEALIZADO
    if option[3] == "Y":
        flow_fline = fline(
            p_gen,
            p_load,
            z_lines,
            z_trx,
            tap_trx,
            bar_type
        )
        bus = [i+1 for i in range(len(flow_fline))]

#---------------------------------------------ESCRITURA DE SALIDA-----------------------------------------
    #Salida de Gauss-Seidel
    if option[0] == "Y":
        sheet = 'RESULTS GS'
        write(np.array(df_bus[:,0]).flatten(),1,3,SALIDA,sheet)
        write(module_vgs,2,3,SALIDA,sheet)
        write(angle_vgs,3,3,SALIDA,sheet)
        write(p_gs,4,3,SALIDA,sheet)
        write(q_gs,5,3,SALIDA,sheet)
        write(p_gen_gs,6,3,SALIDA,sheet)
        write(q_gen_gs,7,3,SALIDA,sheet)
        write(np.array(df_bus[:,6]).flatten(),8,3,SALIDA,sheet)
        write(np.array(df_bus[:,7]).flatten(),9,3,SALIDA,sheet)
        write(iteration,2,1,SALIDA,sheet)
        write(error,4,1,SALIDA,sheet)
    #Flujo por nexos de Gauss-Seidel
        sheet = "POWER FLOW GS"
        write(np.array(df_lines[:,0]).flatten(),1,2,SALIDA,sheet)
        write(np.array(df_trx[:,0]).flatten(),1,len(index_i_lines)+2,SALIDA,sheet)
        write(index_i_lines,2,2,SALIDA,sheet)
        write(index_i_trx,2,len(index_i_lines)+2,SALIDA,sheet)
        write(index_j_lines,3,2,SALIDA,sheet)
        write(index_j_trx,3,len(index_j_lines)+2,SALIDA,sheet)
        write(p_ij,4,2,SALIDA,sheet)
        write(q_ij,5,2,SALIDA,sheet)
        write(p_ji,6,2,SALIDA,sheet)
        write(q_ji,7,2,SALIDA,sheet)
        write(p_losses,8,2,SALIDA,sheet)
        write(q_losses,9,2,SALIDA,sheet)

    #Salida de Newton-Raphson
    if option[1] == "Y":
        sheet = 'RESULTS NR'
        write(np.array(df_bus[:,0]).flatten(),1,3,SALIDA,sheet)
        write(flow_nr_mod,2,3,SALIDA,sheet)
        write(np.degrees(flow_nr_ang),3,3,SALIDA,sheet)
        write(p_nr,4,3,SALIDA,sheet)
        write(q_nr,5,3,SALIDA,sheet)
        write(p_gen_nr,6,3,SALIDA,sheet)
        write(q_gen_nr,7,3,SALIDA,sheet)
        write(np.array(df_bus[:,6]).flatten(),8,3,SALIDA,sheet)
        write(np.array(df_bus[:,7]).flatten(),9,3,SALIDA,sheet)
        write(iteration_nr,2,1,SALIDA,sheet)
        write(error_nr,4,1,SALIDA,sheet)
    #Flujo por nexos de Newton-Raphson
        sheet = "POWER FLOW NR"
        write(np.array(df_lines[:,0]).flatten(),1,2,SALIDA,sheet)
        write(np.array(df_trx[:,0]).flatten(),1,len(index_i_lines)+2,SALIDA,sheet)
        write(index_i_lines,2,2,SALIDA,sheet)
        write(index_i_trx,2,len(index_i_lines)+2,SALIDA,sheet)
        write(index_j_lines,3,2,SALIDA,sheet)
        write(index_j_trx,3,len(index_j_lines)+2,SALIDA,sheet)
        write(p_ij_nr,4,2,SALIDA,sheet)
        write(q_ij_nr,5,2,SALIDA,sheet)
        write(p_ji_nr,6,2,SALIDA,sheet)
        write(q_ji_nr,7,2,SALIDA,sheet)
        write(p_losses_nr,8,2,SALIDA,sheet)
        write(q_losses_nr,9,2,SALIDA,sheet)

    #Salida de Desacoplado Rápido
    if option[2] == "Y":
        sheet = 'RESULTS FD'
        write(np.array(df_bus[:,0]).flatten(),1,3,SALIDA,sheet)
        write(flow_fd_mod,2,3,SALIDA,sheet)
        write(np.degrees(flow_fd_ang),3,3,SALIDA,sheet)
        write(p_fd,4,3,SALIDA,sheet)
        write(q_fd,5,3,SALIDA,sheet)
        write(p_gen_fd,6,3,SALIDA,sheet)
        write(q_gen_fd,7,3,SALIDA,sheet)
        write(np.array(df_bus[:,6]).flatten(),8,3,SALIDA,sheet)
        write(np.array(df_bus[:,7]).flatten(),9,3,SALIDA,sheet)
        write(iteration_fd,2,1,SALIDA,sheet)
        write(error_fd,4,1,SALIDA,sheet)
    #Flujo por nexos de Desacoplado rápido
        sheet = "POWER FLOW FD"
        write(np.array(df_lines[:,0]).flatten(),1,2,SALIDA,sheet)
        write(np.array(df_trx[:,0]).flatten(),1,len(index_i_lines)+2,SALIDA,sheet)
        write(index_i_lines,2,2,SALIDA,sheet)
        write(index_i_trx,2,len(index_i_lines)+2,SALIDA,sheet)
        write(index_j_lines,3,2,SALIDA,sheet)
        write(index_j_trx,3,len(index_j_lines)+2,SALIDA,sheet)
        write(p_ij_fd,4,2,SALIDA,sheet)
        write(q_ij_fd,5,2,SALIDA,sheet)
        write(p_ji_fd,6,2,SALIDA,sheet)
        write(q_ji_fd,7,2,SALIDA,sheet)
        write(p_losses_fd,8,2,SALIDA,sheet)
        write(q_losses_fd,9,2,SALIDA,sheet)

    #Salida de Linealizado
    if option[3] == "Y":
        sheet = 'RESULTS DC'
        write(bus,1,2,SALIDA,sheet)
        write(flow_fline,2,2,SALIDA,sheet)


 

if __name__=="__main__":
    main()

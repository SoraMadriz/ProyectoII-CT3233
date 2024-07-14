import numpy as np      #Dependencia de la estructura del programa
import pandas as pd     #Dependencia de la lectura y escritura

document = "data_io2"
path = "~/Documentos/USB/ProyectoII"

#LECTURA DE HOJA DE CONFIGURACIONES
df_config = np.matrix(pd.read_excel(f'{path}/files/{document}.xlsx',sheet_name="CONFIG"))

#LECTURA DE HOJA DE BUS
df_bus = np.matrix(pd.read_excel(f'{path}/files/{document}.xlsx',sheet_name="BUS"),dtype=object)
#Modelo zip
zip_matrix = df_bus[:,[-1,-2,-3]]
#Ordenar la matrix bus
mask = np.array(df_bus[:, 0] != 'OFF').flatten() #Identifica los buses apagados
df_bus = df_bus[mask]                            #Elimina los buses apagados
df_bus = df_bus[:,1:]                            #Eliminando columnas de STATUS 
buses = df_bus[:,1].astype(int)                  #Buses de conexión como enteros
df_bus[:,1] = buses                              #Conversión de buses a enteros
df_bus = np.delete(df_bus,(-1,-2,-3),axis=1)     #Descargtando valores a no utilizar
#LLenando valores por los nominales(Voltajes mod)
location = np.where(df_bus[:,3] == '-')[0]       #Identificar valores desconocidos
df_bus[location,3] = 1                           #Asignarle el valor de 1
#LLenando valores por los nominales(Voltajes ang)
location = np.where(df_bus[:,4] == '-')[0]       #Identificar valores desconocidos
df_bus[location,4] = 0                           #Asignarle el valor de 0
#LLenando valores por los nominales(P generada)
location = np.where(df_bus[:,5] == '-')[0]       #Identificar valores desconocidos
df_bus[location,5] = 0                           #Asignarle el valor de 0
#LLenando valores por los nominales(Q generada)
location = np.where(df_bus[:,6] == '-')[0]       #Identificar valores desconocidos
df_bus[location,6] = 0                           #Asignarle el valor de 0
#Fasores de voltajes
ang = np.array(df_bus[:,4]).flatten()            #Angulo del voltaje
mod = np.array(df_bus[:,3]).flatten()            #Modulo del voltaje
ang = np.pi*(ang)/180                            #De grados a radianes
df_bus[:,3] = [[i_mod*(np.cos(i_ang)+1j*np.sin(i_ang))] for i_mod,i_ang in zip(mod,ang)]
df_bus = np.delete(df_bus,4,axis=1)
#Nueva potencia según el modelo zip
for locale,values_zip in enumerate(zip_matrix):
    model_zip = np.array(values_zip).flatten()
    if model_zip[0]=="-%":
        continue
    else:
        df_bus[locale,6] = df_bus[locale,6]*(model_zip[0]*(abs(df_bus[locale,3])**2) + model_zip[1]*abs(df_bus[locale,3]) + model_zip[2])
        df_bus[locale,7] = df_bus[locale,7]*(model_zip[0]*(abs(df_bus[locale,3])**2) + model_zip[1]*abs(df_bus[locale,3]) + model_zip[2])



#LECTURA DE HOJA DE LINES
df_lines = np.matrix(pd.read_excel(f'{path}/files/{document}.xlsx',sheet_name="LINES"),dtype=object)
#Ordenar la matrix lines
mask = np.array(df_lines[:, 0] != 'OFF').flatten()   #Identifica las líneas apagadas
df_lines = df_lines[mask]                            #Elimina las líneas apagadas
df_lines = df_lines[:,1:]                            #Eliminando columnas de STATUS
buses = df_lines[:,1:2].astype(int)                  #Buses de conexión como enteros
df_lines[:,1:2] = buses                              #Conversión de buses a enteros
#LLenando valores por los nominales(B SHUNT)
location = np.where(df_lines[:,-1] == '-')[0]        #Identificar valores desconocidos
df_lines[location,-1] = 0                            #Asignarle el valor de 0


#LECTURA DE HOJA DE TRX
df_trx = np.matrix(pd.read_excel(f'{path}/files/{document}.xlsx',sheet_name="TRX"),dtype=object)
#Ordenar la matrix lines
mask = np.array(df_trx[:, 0] != 'OFF').flatten()     #Identifica las TRX apagadas
df_trx = df_trx[mask]                                #Elimina las TRX apagadas
df_trx = df_trx[:,1:]                                #Eliminando columnas de STATUS
buses = df_trx[:,1:3].astype(int)                    #Buses de conexión como enteros
df_trx[:,1:3] = buses                                #Conversión de buses a enteros

#LECTURA DE HOJA DE SHUNT_ELEMENTS
df_sht = np.matrix(pd.read_excel(f'{path}/files/{document}.xlsx',sheet_name="SHUNT_ELEMENTS"),dtype=object)
#Ordenar la matrix SHUNT_ELEMENTS
mask = np.array(df_sht[:, 0] != 'OFF').flatten()     #Identifica las SHUNT ELEMENTS apagadas
df_sht = df_sht[mask]                                #Elimina las SHUNT ELEMENTS apagadas
df_sht = df_sht[:,1:]                                #Eliminando columnas de STATUS
buses = df_sht[:,1].astype(int)                      #Buses de conexión como enteros
df_sht[:,1] = buses                                  #Conversión de buses a enteros

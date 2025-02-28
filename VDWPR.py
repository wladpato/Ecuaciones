
import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import plotly.express as px

from PyInstaller.utils.hooks import copy_metadata

from scipy.optimize import newton
from math import exp

# Configurar la página para usar todo el ancho
st.set_page_config(layout="wide")



datas = copy_metadata("streamlit")


st.subheader("EOS Ecuaciones de Estado y Factor de Compresibilidad Z")
t0, t1, t2, t3, t4, t5 = st.tabs(["Gas Ideal","Isotermas Van der Waals", "Isotermas Peng Robinson", "Factor de Compresibilidad Z", "Dranchuk y Abou-Kassem", "Fundamento Teórico"])
    

# Carga de datos
archivo_excel = 'Datos_PcTc.xlsx'
df1 = pd.read_excel(archivo_excel, sheet_name='Hoja1')
df = df1.copy()


df2 = pd.read_excel(archivo_excel, sheet_name='Hoja2')

# Selección de componente
Components = df['Formula_Química'].unique()

st.sidebar.markdown("## Autor: Wladimir Chávez")
selected_component = st.sidebar.selectbox("**Seleccionar un componente:**", Components)

# Filtrado de datos
df_component = df[df['Formula_Química'] == selected_component]
df_component.reset_index(drop=True, inplace=True)

# Conversión de unidades de presión de psia a Pa
df_component.loc[:,'Pc'] = df_component['Pc'] / 14.7 * 100000

# Constante universal de los gases (R) en (m3.Pa)/(mol.K)
R = 8.314

# Valores críticos del componente seleccionado
Pc = df_component.loc[0, 'Pc']  # Asegúrate de obtener un valor escalar
Tc = df_component.loc[0, 'Tc'] # Asegúrate de obtener un valor escalar
w = df_component.loc[0, 'w'] # Asegúrate de obtener un valor escalar


df_component2 = df2[df2['Formula_Química'] == selected_component]
df_component2.reset_index(drop=True, inplace=True)


######## Temperature Range ########
T1, T2 = -50, 120                         # Start and end temperatures, °C

st. sidebar.write("Rango de temperaturas:")
T1 = st.sidebar.slider("T °C:",-273,120,value=-50)
T2 = st.sidebar.slider("T °C:",-250,120,value=120)

T_step = 10                               # Step size, °C
T = np.arange(T1 + 273.15, T2 + 273.15, T_step)  # Discretisation and temperature conversion en Kelvin

T = np.append(T, Tc)  # para graficar tambien la isoterma critica 


# Obtener listas para Pre, Tem y Vol

Tem_list = df_component2['Tem'].tolist() # solo necesito una temperatura para hacer la isoterma no una lista de Tem   

if len(Tem_list) == 0:
    Tem_list = [0]
else: 
    Tem_Kelv = (Tem_list[0] - 32) * 5/9 + 273.15
    T = np.append(T, Tem_Kelv)  # para graficar tambien la isoterma critica 
Tem = Tem_list[0]
Vol = df_component2['Vol'].tolist()
Pre = df_component2['Pre'].tolist()


######## Molar Volume Range ########
V1, V2 = 0.00006, 0.001     # Start and end molar volume, m3

st. sidebar.write("Rango de Volumen:")
V1 = st.sidebar.slider("V min:",1,100,value=4)
V2 = st.sidebar.slider("Vmax:",6,100,value=100)

V1 = V1/100000
V2 = V2/100000


V_step = 0.000001           # Step size, m3
V = np.arange(V1, V2, V_step)  # Discretisation

######## Pressure Calculation Van der  Waals ########
def vdw(T, V):
    # Substance-specific constants
    a = (27 * R**2 * Tc**2) / (64 * Pc)
    b = (R * Tc) / (8 * Pc)
    P = np.zeros((len(T), len(V)))
    for i in range(len(T)):
        for j in range(len(V)):
            P[i, j] = ((R * T[i]) / (V[j] - b) - a / V[j]**2) / 100000
                   
    return P



P_vdw = vdw(T, V)

# Convierte la presion de Bar a  PSI
P_vdwP= P_vdw * 14.7
# Convierte la temperatura de K a °F
T_F = (T - 273.15) * 9/5 + 32



#  T en kelvin
# Constante universal de los gases (R) en (m3.Pa)/(mol.K) R = 8.314
# Función para calcular el factor de compresibilidad Z
def calcular_Z(T, V, Pc, Tc):
    # Constantes específicas de la sustancia
    a = (27 * R**2 * Tc**2) / (64 * Pc)
    b = (R * Tc) / (8 * Pc)
    
    ##V1, V2 = 0.00006, 0.001     # Start and end molar volume, m3
    V_step = 0.000001  
    
    # Step size, m3
    V = np.arange(V1, V2, V_step)  # Discretisation
    
    
    Z = np.zeros((len(T), len(V)))
    P2 = np.zeros((len(T), len(V)))
    PTZ = np.zeros((len(T) * len(V), 4))

    c = 0  # Inicialización de c
    for i in range(len(T)):
        for j in range(len(V)):
            P2[i, j] = ((R * T[i]) / (V[j] - b) - a / V[j]**2)
            Z[i, j] = P2[i, j] * V[j] / (R * T[i])
            c += 1
            PTZ[c-1, 0] = c
            PTZ[c-1,1]=P2[i, j]
            PTZ[c-1,2]=T[i]
            PTZ[c-1,3]=Z[i, j]
    return PTZ





PTZ = calcular_Z(T, V, Pc, Tc)


with t0:


    def ley_gases_ideales(T, V):
        P_ideal = np.zeros((len(T), len(V)))
        for i in range(len(T)):
            for j in range(len(V)):
                P_ideal[i, j] = (R * T[i]) / V[j]*14.7/ 100000
        return P_ideal


    P_idea = ley_gases_ideales(T, V)


    #st.subheader("Gas Ideal:")
    
    # Crear la figura con Plotly
    fig = go.Figure()
    
    
    
    
    
    
    
    Pre = np.zeros((1, len(df2)), dtype=object)
    Pre[0, :] = df2.iloc[:, 0]

    
    for i, yi in enumerate(Pre):
        fig.add_trace(go.Scatter(x=Vol, y=yi, mode='lines', name=f'{Tem:.1f} °F(Real)'))

    # Actualizar el layout de la figura para añadir títulos y etiquetas
    fig.update_layout(
        title="Isotermas de Gas Ideal",
        xaxis_title="Volumen (L/mol)",
        yaxis_title="Presión (Psi)",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1,range=[0, 1]),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1000,range=[-1000, 8000]),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad

    ) 
    
    
    
    
    
    
    
    
    for i, yi in enumerate(P_idea):
        fig.add_trace(go.Scatter(x=V * 1000, y=yi, mode='lines', name=f'{T_F[i]:.1f} °F Tc'))

    # Actualizar el layout de la figura para añadir títulos y etiquetas
    fig.update_layout(
        #title="Isotermas de la EOS de Van der Waals",
        xaxis_title="Volumen (L/mol)",
        yaxis_title="Presión (Psi)",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1,range=[0, 1]),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1000,range=[-1000, 8000]),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad

    ) 

    # Utiliza Streamlit para mostrar la figura interactiva
    st.plotly_chart(fig, use_container_width=True)





with t1:

    
    Pre = np.zeros((1, len(df2)), dtype=object)
    Pre[0, :] = df2.iloc[:, 0]

    # Grafico de la isoterma real 
    fig = go.Figure()
    
    for i, yi in enumerate(Pre):
        fig.add_trace(go.Scatter(x=Vol, y=yi, mode='lines', name=f'{Tem:.1f} °F(Real)'))

    # Actualizar el layout de la figura para añadir títulos y etiquetas
    fig.update_layout(
        title="Isotermas de Van der Waals",
        xaxis_title="Volumen (L/mol)",
        yaxis_title="Presión (Psi)",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1,range=[0, 1]),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1000,range=[-1000, 8000]),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad

    ) 
        

    for i, yi in enumerate(P_vdwP):
        fig.add_trace(go.Scatter(x=V * 1000, y=yi, mode='lines', name=f'{T_F[i]:.1f} °F'))

    # Actualizar el layout de la figura para añadir títulos y etiquetas
    fig.update_layout(
        #title="Isotermas de la EOS de Van der Waals",
        xaxis_title="Volumen (L/mol)",
        yaxis_title="Presión (Psi)",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1,range=[0, 1]),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1000,range=[-1000, 8000]),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad

    ) 

    # Utiliza Streamlit para mostrar la figura interactiva
    #st.plotly_chart(fig, use_container_width=True)
   
   
       # Crear la figura con Plotly    datos exp
    #fig = go.Figure()


    # Utiliza Streamlit para mostrar la figura interactiva
    st.plotly_chart(fig, use_container_width=True)
   
    
    # Visualización de la tabla con los valores de Pc y Tc con unidades
    df_component.loc[:,'Pc'] = df_component['Pc'] * 14.7 / 100000
    #df_component['Tc'] = (df_component['Tc'] − 273.15) × 9/5 + 32
    df_component.loc[:,'Pc_'] = df_component['Pc'].apply(lambda x: f"{x:.2f} PSI")
    df_component.loc[:,'Tc_'] = df_component['Tc'].apply(lambda x: f"{x:.2f} °K")
    df_component.loc[:,'Tc_F'] = (df_component['Tc'] - 273.15) * 9/5 + 32
    df_component.loc[:,'Tc_F'] = df_component['Tc_F'].apply(lambda x: f"{x:.2f} °F")
    
    b = (R * Tc) / (8 * Pc)*1000
    df_component.loc[:,'b Covolumen'] = b
    
    
    st.write("Valores Críticos del Componente Seleccionado:")
    st.table(df_component[['Formula_Química', 'Pc_', 'Tc_','Tc_F','b Covolumen']])
    
    #Grafico 3d
    #X, Y = np.meshgrid(V*1000, T)
    #Z = P_vdw


    #fig = go.Figure(data=[go.Surface(z=Z, x=Y, y=X, colorscale='Viridis')])

    # Actualizar el diseño de la figura
    #fig.update_layout(
        #title="Ecuación de Estado de Van der Waals",
        #font=dict(color="Black"),
        #width=800,
        #height=600,
        #margin=dict(l=0, r=0, b=0, t=50),
        #scene=dict(
            #yaxis=dict(nticks=10, range=[0, 0.6]),
            #xaxis=dict(nticks=T_step, range=[T.max(), T.min() + T_step]),
            #zaxis=dict(nticks=10, range=[P_vdw.min(), P_vdw.max()]),
            #yaxis_title="Volumen (L) Eje:Y",
            #xaxis_title='Temperatura (K) Eje:X',
            #zaxis_title='Presión (bar) Eje:Z',
            #bgcolor='white'  # Cambia el fondo a blanco
        #)
    #)

    # Utiliza Streamlit para mostrar la figura interactiva
    #st.plotly_chart(fig, use_container_width=True)
    
    
    
    
    
    

with t2:


    ######## Pressure Calculation Peng Robinson########
    def vdw(T, V):
        # Substance-specific constants
        a = ( 0.45724* R**2 * Tc**2) / (Pc)
        b = (0.07780*R * Tc) / (Pc)
        P = np.zeros((len(T), len(V)))
    
        for i in range(len(T)):
            for j in range(len(V)):
                alpha = (1 + (0.37464 + 1.54226*w - 0.26992*w**2)*(1-(T[i]/Tc)**0.5))**2
                P[i, j] = ((R * T[i]) / (V[j] - b) - alpha*a /((V[j]**2 + 2*V[j]*b - b**2)) ) / 100000
                   
        return P
    P_vdw = vdw(T, V)

    # Convierte la presion de Bar a  PSI
    P_vdwP= P_vdw * 14.7
    # Convierte la temperatura de K a °F
    T_F = (T - 273.15) * 9/5 + 32
    
       # Crear la figura con Plotly
    fig = go.Figure()
    
    
        # Crear la figura datos reales
    fig = go.Figure()
    #yi = np.array([yi])
    for i, yi in enumerate(Pre):
        fig.add_trace(go.Scatter(x=Vol, y=yi, mode='lines', name=f'{Tem:.1f} °F(Real)'))
    # Actualizar el layout de la figura para añadir títulos y etiquetas
    fig.update_layout(
        title="Isotermas de Peng Robinson",
        xaxis_title="Volumen (L/mol)",
        yaxis_title="Presión (Psi)",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1,range=[0, 1]),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1000,range=[-1000, 8000]),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad

    ) 
    
    
    
    

    for i, yi in enumerate(P_vdwP):
        fig.add_trace(go.Scatter(x=V * 1000, y=yi, mode='lines', name=f'{T_F[i]:.1f} °F'))

    # Actualizar el layout de la figura para añadir títulos y etiquetas
    fig.update_layout(
        #title="Isotermas de la EOS de Van der Waals",
        xaxis_title="Volumen (L/mol)",
        yaxis_title="Presión (Psi)",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1,range=[0, 1]),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1000,range=[-1000, 8000]),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad

    ) 

    # Utiliza Streamlit para mostrar la figura interactiva
    st.plotly_chart(fig, use_container_width=True)
    
    # Crear la figura con Plotly
    fig = go.Figure()
       
    
    # Visualización de la tabla con los valores de Pc y Tc con unidades
    df_component.loc[:,'Pc'] = df_component['Pc'] 
    #df_component['Tc'] = (df_component['Tc'] − 273.15) × 9/5 + 32
    df_component.loc[:,'Pc_'] = df_component['Pc'].apply(lambda x: f"{x:.2f} PSI")
    df_component.loc[:,'Tc_'] = df_component['Tc'].apply(lambda x: f"{x:.2f} °K")
    
    df_component.loc[:,'Tc_F'] = (df_component['Tc'] - 273.15) * 9/5 + 32
    df_component.loc[:,'Tc_F'] = df_component['Tc_F'].apply(lambda x: f"{x:.2f} °F")

  
    
    df_component.loc[:,'w_Factor acéntrico'] = df_component['w'].apply(lambda x: f"{x:.3f} ")
    
    
    
    df_component.loc[:,'b Covolumen'] = b
    
    
    st.write("Valores Críticos del Componente Seleccionado:")
    st.table(df_component[['Formula_Química', 'Pc_', 'Tc_','Tc_F','w_Factor acéntrico','b Covolumen'] ])
    
    #Grafico 3d
    #X, Y = np.meshgrid(V*1000, T)
    #Z = P_vdw


    #fig = go.Figure(data=[go.Surface(z=Z, x=Y, y=X, colorscale='Viridis')])

    # Actualizar el diseño de la figura
    #fig.update_layout(
        #title="Ecuación de Estado de Van der Waals",
        #font=dict(color="Black"),
        #width=800,
        #height=600,
        #margin=dict(l=0, r=0, b=0, t=50),
        #scene=dict(
            #yaxis=dict(nticks=10, range=[0, 0.6]),
            #xaxis=dict(nticks=T_step, range=[T.max(), T.min() + T_step]),
            #zaxis=dict(nticks=10, range=[P_vdw.min(), P_vdw.max()]),
            #yaxis_title="Volumen (L) Eje:Y",
            #xaxis_title='Temperatura (K) Eje:X',
            #zaxis_title='Presión (bar) Eje:Z',
            #bgcolor='white'  # Cambia el fondo a blanco
        #)
    #)

    # Utiliza Streamlit para mostrar la figura interactiva
    #st.plotly_chart(fig, use_container_width=True)
    
    
    
    




with t5:
    st.markdown("""
        # Las Isotermas de Van der Waals
        Son curvas que representan cómo se comportan los gases reales en función de la presión y el volumen a una temperatura constante. Estas isotermas se derivan de la ecuación de van der Waals, que es una modificación de la ley de los gases ideales para tener en cuenta el tamaño no nulo de las partículas y las fuerzas de atracción e interacciones entre ellas.
        
        ## Ecuación:
        $$
        \\left( P + \\frac{{n^2 \\cdot a}}{{V^2}} \\right) \\left( V - n \\cdot b \\right) = n \\cdot R \\cdot T
        $$
        
        Donde:
        - \(P\) es la presión del gas.
        - \(V\) es el volumen.
        - \(n\) es el número de moles.
        - \(R\) es la constante de los gases ideales.
        - \(T\) es la temperatura absoluta.
        - \(a\) y \(b\) son constantes características de cada gas.
        
        El término “a” toma en cuenta las fuerzas atractivas entre partículas de Van der Waals
        
        El término “b” toma en cuenta el volumen finito de las moléculas es decir que las moléculas ocupan un volumen real. 
        
        $$       
        a = \\frac{27 R^2 Tc^2}{64 Pc}
        $$ 
        
        $$ 
        b = \\frac{R Tc}{8 Pc}
        $$

        
        Las isotermas de un gas real tienen una forma más compleja que las isotermas de un gas ideal (hipérbolas), ya que deben dar cuenta de los cambios de fase que puede experimentar.

        En la figura inferior se han representado las denominadas isotermas de Andrews. Dichas isotermas fueron medidas experimentalmente, y representan la presión en función del volumen a distintas temperaturas
        
        La isoterma representada en rojo se denomina isoterma crítica (y su temperatura, la temperatura crítica). Esta isoterma separa dos comportamientos: cuando una sustancia se encuentra a una temperatura superior a su temperatura crítica, siempre está en estado gaseoso, por muy alta que sea la presión. Por el contrario, cuando está a una temperatura inferior a la crítica, puede estar en estado sólido, líquido o vapor (en la gráfica se han representado solamente las zonas de líquido y vapor).
        """)
    
    from PIL import Image
    # Carga la imagen desde el archivo "imagen.png"
    imagen = Image.open("Isotermas.png")

    # Muestra la imagen en Streamlit
    st.image(imagen, caption="Imagen de ejemplo", use_column_width=True, width=500)
    
    st.markdown("""
        # Factor de compresibilidad
        Es un coeficiente que describe la relación entre el volumen molar de un gas real y el volumen molar de un gas ideal a la misma temperatura y presión.
        
        $$       
        Z = \\frac{P T}{n R T}
        $$ 
        
    """)



# Agregar una nueva pestaña para visualizar el factor Z
#with st.tabs("Factor Z"):
with t3:   
    
    R = 8.314  # Gas constant

    def calcular_Z(T, V, Pc, Tc):
        # Constantes específicas de la sustancia
        a = (27 * R**2 * Tc**2) / (64 * Pc)
        b = (R * Tc) / (8 * Pc)
    
        ##V1, V2 = 0.00006, 0.001     # Start and end molar volume, m3
        V_step = 0.000001            # Step size, m3
        V = np.arange(V1, V2, V_step)  # Discretisation
    
        Z = np.zeros((len(T), len(V)))
        P2 = np.zeros((len(T), len(V)))
        P2psi = np.zeros((len(T), len(V)))
    
        for i in range(len(T)):
            for j in range(len(V)):
                P2[i, j] = ((R * T[i]) / (V[j] - b) - a / V[j]**2)
                P2psi[i, j] = ((R * T[i]) / (V[j] - b) - a / V[j]**2) *14.7/ 100000
                Z[i, j] = P2[i, j] * V[j] / (R * T[i])
    
        return P2psi, T, Z




    T = np.linspace(300, 400, 10)  # Temperature range

    P2psi, T_grid, Z_grid = calcular_Z(T, V, Pc, Tc)

    # Create the surface plot
    #surface = go.Surface(z=Z_grid, x=T_grid, y=P2psi)
    surface = go.Surface(z=Z_grid, x=P2psi, y=T_grid)
    data = [surface]

    layout = go.Layout(
        title='Factor de compresibilidad (Z)',
        font=dict(color="Black"),
        width=800,
        height=600,
        margin=dict(l=0, r=0, b=0, t=50),
        scene=dict(
            xaxis=dict(title='Presion (psia) Eje:X'),
            yaxis=dict(title='Temperatura (K) Eje:Y'),
            zaxis=dict(title='Factor de compresibilidad (Z) Eje:Z'),
            bgcolor='white'  # Cambia el fondo a blanco
        )
    )
    # Configura el eje X para mostrar los valores sin prefijos y sin redondear
    layout.scene.xaxis.tickformat = '.0f'  # Dos decimales
    fig = go.Figure(data=data, layout=layout)

    # Display the figure with Streamlit
    st.plotly_chart(fig)
        
    from mpl_toolkits.mplot3d import Axes3D
    import plotly.express as px

    # Simulated data (reemplaza esto con tus datos PTZ reales)
    # PTZ = np.random.rand(100, 4)  # Ejemplo: 100 filas y 4 columnas

    # Extract relevant columns
    P = PTZ[:, 1]*14.7/100000
    T = PTZ[:, 2]
    Z = PTZ[:, 3]

    # Create a 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(T, P, Z, c='b', marker='o', label='PTZ')

    # Customize plot labels
    ax.set_xlabel('Temperatura (K)')
    ax.set_ylabel('Presión (psia)')
    ax.set_zlabel('Compressibility Factor (Z)')

    # Show the plot
    st.pyplot(fig)  # Muestra el gráfico en Streamlit
   
   
   



with t4:

    # Definimos los coeficientes como constantes Dranchuk y Abou-Kassem
    A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11 = 0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210
    
    # Definir coeficientes dependientes de Tr
    def coefficients(Tr):
        C1 = A1 + A2 / Tr + A3 / Tr**3 + A4 / Tr**4 + A5 / Tr**5
        C2 = A6 + A7 / Tr + A8 / Tr**2
        C3 = (A7 / Tr + A8 / Tr**2)
        C4 = A10 * (1 + A11 / Tr**2) / Tr**3
        return C1, C2, C3, C4
    
    # Función f(Z) que depende de Z, Pr, Tr y ρr
    def f_z(Z, Pr, Tr):
        C1, C2, C3, C4 = coefficients(Tr)
        rho_r = 0.27 * Pr / (Z * Tr)
        fZ = Z - (1 + C1 * rho_r + C2 * rho_r**2 -  A9 * C3 * rho_r**5 + A10 * (1 + A11 * rho_r**2) * (rho_r**2 / Tr**3) * exp(-A11 * rho_r**2))
        return fZ
    
    # Método de Newton-Raphson utilizando scipy.optimize.newton para encontrar Z
    def calculate_z(Pr, Tr, tol=1e-6, max_iter=500, smart_guess=True, newton_kwargs=None):
        guess = 0.9 if Pr < 15 else 2.0  # Adivinanza inicial básica
        
        def _check_working_Pr_Tr_range(Pr, Tr):
            return Pr < 15 and 1.0 < Tr < 3.0  # Ejemplo de rango de operación
        
        def _get_z_model(model='kareem'):
            if model == 'kareem':
                return lambda Pr, Tr: 1.0 - (0.01 * Pr)  # Ejemplo simplificado de un modelo explícito
        
        worked = False
        guesses = []
        
        if smart_guess:
            if _check_working_Pr_Tr_range(Pr, Tr):
                smart_guess_model = 'kareem'
                guess_zmodel_func = _get_z_model(model=smart_guess_model)
                guess_smart = guess_zmodel_func(Pr=Pr, Tr=Tr)
                guesses.append(guess_smart)
            guesses.append(guess)
        else:
            guesses.append(guess)
        
        # Adivinanzas adicionales si la primera no funciona
        guesses.extend([guess * 0.5, guess * 1.5, 1.0])
        
        for guess_ in guesses:
            try:
                if newton_kwargs is None:
                    Z_solution = newton(f_z, guess_, args=(Pr, Tr), tol=tol, maxiter=max_iter)
                else:
                    Z_solution = newton(f_z, guess_, args=(Pr, Tr), tol=tol, maxiter=max_iter, **newton_kwargs)
                worked = True
            except RuntimeError:
                continue
            
            if worked:
                if abs(f_z(Z_solution, Pr, Tr)) < tol:
                    return Z_solution
        
        raise RuntimeError("No se alcanzó la convergencia después de probar múltiples conjeturas.")
    
    # Interfaz de usuario
    st.title("Gráfico del Factor de Compresibilidad Z (Dranchuk y Abou-Kassem)")
    Pr_range = st.slider("Selecciona el rango de Presión Pseudo-reducida (Pr)", 0.05, 15.0, (0.1, 8.0))
    Tr_values = [1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    
    # Calcular valores de Z para cada Tr seleccionado
    Pr_values = np.linspace(Pr_range[0], Pr_range[1], 100)
    Z_results = {}
    
    for Tr in Tr_values:
        Z_values = [calculate_z(Pr, Tr) for Pr in Pr_values]
        Z_results[Tr] = Z_values
    
    # Crear el gráfico interactivo con Plotly
    fig = go.Figure()

    for Tr, Z_values in Z_results.items():
        fig.add_trace(go.Scatter(x=Pr_values, y=Z_values, mode='lines', name=f'Tr = {Tr:.2f}'))
    fig.update_layout(
        title="Factor de Compresibilidad Z vs. Presión Pseudo-reducida Pr",
        xaxis_title="Presión Pseudo-reducida Pr",
        yaxis_title="Factor de Compresibilidad Z",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad
    )

    st.plotly_chart(fig)

    ######## CALCULO DEL VOLUMEN A PARTIR DE A EC.PV=znRT OJO: al usar Dranchuk y Abou-Kassem para calcular Z solo se puede usar  para  Temperaturas mayores a la Tc ya que Tr=T/Tc debe ser mayor a 1########

    Pr2 = np.arange(0.4, 12, 0.045) # OJO pasos muy pequeños pueden afectar la convergencia despues para encontrar Z 


    #Ts = np.arange(310, 340, 10)  # grados K
    Ts = np.arange(Tc + 5.79, Tc+100, 10)  # grados K
    Ts = np.append(Ts, 373.15)
    Ts_F = (Ts - 273.15) * 9/5 + 32  # grados F

    Tr2 = Ts / Tc   # K/K
    
    Pss = Pr2 * Pc  # Pascal por que Pc  esta convertido en pascal, ajusta según corresponda
    P2_pascal = Pss 
    Pss_psi=Pss * 14.7 / 100000


    
    def V_znRT(Ts, P2_pascal):
        V = np.zeros((len(Ts), len(Pr2)))  # Matriz para almacenar volúmenes para cada combinación de Ts y Pr
        for ii in range(len(Ts)):
            for j in range(len(Pr2)):
                Tr = Ts[ii] / Tc  # Relación Tr para cada temperatura
                if 1.0 < Tr < 3.0:  # Condición para el cálculo de Z
                    Z = calculate_z(Pr=Pr2[j], Tr=Tr)
                else:
                    Z = 0  # Si Tr está fuera del rango, se asigna Z=1 para evitar errores
                
                V[ii, j] = Z * 8.314 * Ts[ii] * 1000 / P2_pascal[j]  # Litros, R=8.314 J/mol.K (ejemplo)
        return V
    

    # Como `T` ahora es un array, necesitamos asegurarnos de que el cálculo de V se haga para cada temperatura
    #V_litros = np.array([V_znRT(Ts[i], P2_pascal, Z) for i in range(len(Ts))])
    #V_litros =V_znRT(Ts, P2_pascal, Z)
    V_litros =V_znRT(Ts, P2_pascal)
    
    


    # Crear el gráfico interactivo con Plotly
    fig = go.Figure()

    for i in range(len(Ts)):
        fig.add_trace(go.Scatter(x=V_litros[i], y=Pss_psi, mode='lines', name=f'{Ts_F[i]:.1f} °F'))

    fig.update_layout(
        title="Volumen V=nRTZ/P,  Z mediante método de Dranchuk y Abou-Kassem",
        xaxis_title="Volumen (L/mol)",
        yaxis_title="Presión (Psi)",
        hovermode="closest",
        height=810,
        width=1500,
        xaxis=dict(showgrid=True, gridcolor='lightgray', dtick=0.1,range=[0, 1]),  # Cuadrícula en el eje x
        yaxis=dict(showgrid=True, gridcolor='lightgray', dtick=1000,range=[0, 8000]),  # Cuadrícula en el eje y
        plot_bgcolor='white'  # Fondo blanco para mayor visibilidad

    )   

    #st.plotly_chart(fig)    
    st.plotly_chart(fig, use_container_width=True)
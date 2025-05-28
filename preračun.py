# -*- coding: utf-8 -*-
"""
Created on Sat May 24 00:47:01 2025

@author: Jernej Kušar
"""


#importing paths for lybraries
import sys
sys.path.append("C:/Users/Jernej Kusar/Documents/LFDT splošno/Dodiplomska/Main")


#common libraries
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 1000
import numpy as np
from CoolProp.CoolProp import PropsSI
from pyfluids import Fluid, FluidsList


np.__version__

#custom libraries
import datalib as jdl
import plotlib as jpl



#Uvoz podatkov temperatur (ARSO arhiv)

raw_data = jdl.get_data_fromTXT(r"C:\Users\Jernej Kusar\Documents\ŠOLA 2024-2025 2. semester\Hladilna tehnika in toplotne črpalke\Lab. vaje\seminarska_2\dnevna_temp.txt")

raw_data = raw_data[2:]


temp = {}
day = {}
j = 0
for i in range(len(raw_data)):
    if raw_data[i] != "":
        j += 1

        t = float(raw_data[i].split(",")[-1])
        temp.update({f"day-{j}" : t})
        day.update({f"day-{j}" : j})


names = ["day"]


color_palett = {"day" : "red"}
symbol_palett = {"day" : "."}



jpl.plot_xy_data(nozzles=names, values_x=day, values_y=temp, 
                 x_label=r"Dan", y_label=r"Zunanja temperatura [°C]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50,
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = ["day"], connect_line_width = 1,
                 minor_ticks=True, color_palett = color_palett, symbol_palett = symbol_palett, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = True,
                 curve_points = [[[1, 365], [12, 12]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")


#Splošni podatki
T_i_pr = 22 #°C
T_i = 22 #°C
T_ogr_vst_pr = 55 #°C
T_ogr_izt_pr = 45 #°C
T_e_pr = -13 #°C
m = 1.3
izkoristek = 0.75
T_out_pr = -13 #°C
P_v = 50 * 10**-3 #kW
p_okolice = 101325 #Pa
UA = 213 #W/K
T_omejitev = 12 #°C


#Pregretja
D_t_kond = 5 #K
D_t_up = 10 #K
D_t_ogrevanje = 5 #K

#Zrak-voda
D_t_uparjanje_zrak = 7 #K

#Voda-voda
T_podtalnice = 9 #°C
D_t_uparjanje_voda = 3 #K
izkoristek_črpalke = 0.6 #/
D_cevi = 25 * 10**-3 #m
epsilon = 5 * 10**-5 #/
L_cevi = 20 #m
f_D = 0.024336 #https://www.omnicalculator.com/physics/friction-factor #Re = 76288.1, k = 0.00125 mm, D = 25 mm
g = 9.81 #m/s^2

#Izris grafov
krivulja_dan = 50



#Potrebna temperatura na vstopu v ogrevalni sistem

def temp_ogr_vst(T_i_pr, T_i, T_ogr_vst_pr, T_ogr_izt_pr, T_e_pr, T_e, m):
    a = (T_ogr_vst_pr - T_ogr_izt_pr)/2
    b = (T_i_pr - T_e)/(T_i_pr - T_e_pr)
    c = ((T_ogr_vst_pr + T_ogr_izt_pr/2) - T_i)
    d = (T_i_pr - T_e)/(T_i_pr - T_e_pr)

    
    T_ogr_vst = T_i_pr + a * b + c * d**(1/m)

    return T_ogr_vst

vstopna_temp = {}
for i in list(temp.keys()):
    T_e = float(temp[i])
    T_ogr = temp_ogr_vst(T_i_pr, T_i, T_ogr_vst_pr, T_ogr_izt_pr, T_e_pr, T_e, m)
    vstopna_temp.update({i : T_ogr})



mejna_temp = temp_ogr_vst(T_i_pr, T_i, T_ogr_vst_pr, T_ogr_izt_pr, T_e_pr, T_omejitev, m)

jpl.plot_xy_data(nozzles=names, values_x=day, values_y=vstopna_temp, 
                 x_label=r"Dan", y_label=r"Temperatura na vstopu v ogrevalni sistem [°C]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = ["day"], connect_line_width = 1,
                 minor_ticks=True, color_palett = color_palett, symbol_palett = symbol_palett, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = True,
                 curve_points = [[[1, 365], [mejna_temp, mejna_temp]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")



#Izračun podatkov toplotne črpalke zrak-voda

temp_uparjalnika_zrak = np.array(list(temp.values())) + D_t_uparjanje_zrak
temp_uparjalnika_voda = T_podtalnice + D_t_uparjanje_voda

temp_kondenzatorja = np.real(np.array(list(vstopna_temp.values()))) + D_t_ogrevanje

#izračun tlakov

uparjalnik_tlak_zrak = []
kondenzator_tlak = []
for i in range(len(temp)):
    
    kondenzator = Fluid(FluidsList.R290).dew_point_at_temperature(temp_kondenzatorja[i])
    uparjalnik_z = Fluid(FluidsList.R290).bubble_point_at_temperature(temp_uparjalnika_zrak[i])

    p_k = kondenzator.pressure
    p_u_z = uparjalnik_z.pressure

    uparjalnik_tlak_zrak.append(p_u_z)
    kondenzator_tlak.append(p_k)
    
uparjalnik_v = Fluid(FluidsList.R290).bubble_point_at_temperature(temp_uparjalnika_voda)
uparjalnik_tlak_voda = uparjalnik_v.pressure


#Izri krivulje entalpija-tlak

range_tlak = np.linspace(0.5, 100, 1000) #bar (več kot 0)

krivulja_tlak_1 = {}
krivulja_tlak_2 = {}

krivulja_entalpija_levo = {}
krivulja_entalpija_desno = {}
j = 0
for i in range_tlak:
    try:
        levo = Fluid(FluidsList.R290).dew_point_at_pressure(i*10**5).enthalpy * 10**-3
        desno = Fluid(FluidsList.R290).bubble_point_at_pressure(i*10**5).enthalpy * 10**-3
        j+=1
        krivulja_tlak_1.update({f"R290-{j}" : i})
        krivulja_tlak_2.update({f"R290_1-{j+len(range_tlak)}" : i})
        krivulja_entalpija_levo.update({f"R290-{j}" : levo})
        krivulja_entalpija_desno.update({f"R290_1-{j+len(range_tlak)}" : desno})
    except Exception as e:
        print(f"Tlak vrha: {i}")
        print(f"Exception: {e}")
        break

krivulja_entalpija = krivulja_entalpija_levo | krivulja_entalpija_desno
krivulja_tlak = krivulja_tlak_1 | krivulja_tlak_2


krivulja = ["R290", "R290_1"]
kc_pal = {"R290": "green", "R290_1" : "green"}
ks_pal = {"R290": ".", "R290_1" : "."}



#Izračun, katere dneve ogrevamo
ogrevani_dnevi = []
for i in list(temp.keys()):
    if temp[i] <= 12:
        ogrevani_dnevi.append(i)
    else:
        ogrevani_dnevi.append(None)


#parametri toplotne črpalke zrak-voda


def h_2_real(P_up, P_kond, T_up, ni):
    h1 = PropsSI("H", "P", P_up, "T", T_up, "R290") * 10**-3
    s1 = PropsSI("S", "P", P_up, "T", T_up, "R290")
    h2_id = PropsSI("H", "P", P_kond, "S", s1, "R290") * 10**-3   
    h2_real = ((h2_id - h1)/ni) + h1
    
    return h2_real


h_zrak = {}

for i in range(len(temp)):
    h3 = PropsSI("H", "P", kondenzator_tlak[i], "T", temp_kondenzatorja[i] - D_t_kond + 273, "R290") * 10**-3
    h4 = h3
    h1 = PropsSI("H", "P", uparjalnik_tlak_zrak[i], "T", temp_uparjalnika_zrak[i] + D_t_up + 273, "R290") * 10**-3
    h2_real = h_2_real(uparjalnik_tlak_zrak[i], kondenzator_tlak[i], temp_uparjalnika_zrak[i] + D_t_up + 273, izkoristek)
    
    dan = f"day-{i+1}"
    h_zrak.update({dan : [h1, h2_real, h3, h4]})

h_voda = {}

for i in range(len(temp)):
    h3 = PropsSI("H", "P", kondenzator_tlak[i], "T", temp_kondenzatorja[i] - D_t_kond + 273, "R290") * 10**-3
    h4 = h3
    h1 = PropsSI("H", "P", uparjalnik_tlak_voda, "T", temp_uparjalnika_voda + D_t_up + 273, "R290") * 10**-3
    h2_real = h_2_real(uparjalnik_tlak_voda, kondenzator_tlak[i], temp_uparjalnika_voda + D_t_up + 273, izkoristek)
    
    dan = f"day-{i+1}"
    h_voda.update({dan : [h1, h2_real, h3, h4]})
    

#Izris grafa

temp[f"day-{krivulja_dan}"]

k_x_z = h_zrak[f"day-{krivulja_dan}"]
k_y_z = [uparjalnik_tlak_zrak[list(temp.values()).index(temp[f"day-{krivulja_dan}"])]*10**-5, kondenzator_tlak[list(temp.values()).index(temp[f"day-{krivulja_dan}"])]*10**-5]

c_p_z = [[(k_x_z[-1], k_x_z[0]), (k_y_z[0], k_y_z[0])], [(k_x_z[0], k_x_z[1]), (k_y_z[0], k_y_z[1])], [(k_x_z[1], k_x_z[2]), (k_y_z[1], k_y_z[1])], [(k_x_z[2], k_x_z[-1]), (k_y_z[1], k_y_z[0])]]

jpl.plot_xy_data(nozzles=krivulja, values_y=krivulja_tlak, values_x=krivulja_entalpija, 
                 y_label=r"Tlak [bar]", x_label=r"Entalpija [kJ/kg]", title=f"Krivulja zrak-zrak za dan: {krivulja_dan}", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 marker_size = 50, scale_y = "log", connect_points_sort = "y",
                 show_grid = True, connect_line_color = {"R290": "green", "R290_1":"green"},
                 point_opacity = 0, connect_points = ["R290", "R290_1"], connect_line_width = 1,
                 minor_ticks=True, color_palett = kc_pal, symbol_palett = ks_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = True,
                 curve_points = c_p_z, curve_color = "violet", curve_width = 1, curve_style="-",
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 9,
                 title_fontsize = 9, plot_font = "Palatino Linotype")

k_x = h_voda[f"day-{krivulja_dan}"]
k_y = [uparjalnik_tlak_voda*10**-5, kondenzator_tlak[list(temp.values()).index(temp[f"day-{krivulja_dan}"])]*10**-5]

c_p_v = [[(k_x[-1], k_x[0]), (k_y[0], k_y[0])], [(k_x[0], k_x[1]), (k_y[0], k_y[1])], [(k_x[1], k_x[2]), (k_y[1], k_y[1])], [(k_x[2], k_x[-1]), (k_y[1], k_y[0])]]


jpl.plot_xy_data(nozzles=krivulja, values_y=krivulja_tlak, values_x=krivulja_entalpija, 
                 y_label=r"Tlak [bar]", x_label=r"Entalpija [kJ/kg]", title=f"Krivulja voda-zrak za dan: {krivulja_dan}", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 marker_size = 50, scale_y = "log", connect_points_sort = "y",
                 show_grid = True, connect_line_color = {"R290": "green", "R290_1":"green"},
                 point_opacity = 0, connect_points = ["R290", "R290_1"], connect_line_width = 1,
                 minor_ticks=True, color_palett = kc_pal, symbol_palett = ks_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = True,
                 curve_points = c_p_v, curve_color = "violet", curve_width = 1, curve_style="-",
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 9,
                 title_fontsize = 9, plot_font = "Palatino Linotype")

#Potrebna toplotna moč


def izgube(UA, T_zunaj, T_notranja):
    Q = UA * (T_notranja - T_zunaj)
    return Q


potrebno_ogrevanje = {}
for key, val in temp.items():
    q = izgube(UA, val, T_i)
    potrebno_ogrevanje.update({key: q})


#Masni pretok R290 glede na zunanje pogoje

masni_pretok_voda = {}
masni_pretok_zrak = {}
for i in temp.keys():
    Q_izgube= izgube(UA, temp[i], T_i)
    
    h_v_2, h_v_3 = h_voda[i][1], h_voda[i][2]
    h_z_2, h_z_3 = h_zrak[i][1], h_zrak[i][2]
    
    mf_zrak = - Q_izgube/((h_z_3 - h_z_2)*10**3)
    mf_voda = - Q_izgube/((h_v_3 - h_v_2)*10**3)
    
    masni_pretok_voda.update({i: mf_voda})
    masni_pretok_zrak.update({i : mf_zrak})
    
    
#COP zrak

Q_up_zrak = {}
Q_kond_zrak = {}
P_komp_zrak = {}
COP_zrak = {}
dnevi = {}
j = 0
for i in list(temp.keys()):
    Q_up = masni_pretok_zrak[i] * (h_zrak[i][0] - h_zrak[i][-1])
    Q_kond = - masni_pretok_zrak[i] * (h_zrak[i][2] - h_zrak[i][1])
    P_komp = masni_pretok_zrak[i] * (h_zrak[i][1] - h_zrak[i][0])
    
    COP = Q_kond / (P_komp + P_v)
    j+=1
    
    if i in ogrevani_dnevi:
        dnevi.update({i : j})
        P_komp_zrak.update({i : P_komp})
        Q_up_zrak.update({i : Q_up})
        Q_kond_zrak.update({i : Q_kond})
        COP_zrak.update({i : COP})
        
    else:
        n = f"negreje-{j}"
        dnevi.update({n : j})
        P_komp_zrak.update({n : P_komp})
        Q_up_zrak.update({n : Q_up})
        Q_kond_zrak.update({n : Q_kond})
        COP_zrak.update({n : COP})


denotions = ["day", "negreje"]
c_pal = {"day": "red", "negreje" : "blue"}
s_pal = {"day": ".", "negreje" : "X"}


jpl.plot_xy_data(nozzles=names, values_x=dnevi, values_y=COP_zrak, 
                 x_label=r"Dan", y_label=r"COP [/]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [mejna_temp, mejna_temp]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")



jpl.plot_xy_data(nozzles=denotions, values_x=dnevi, values_y = Q_kond_zrak, 
                 x_label=r"Dan", y_label=r"Toplotni tok kondenzatorja [kW]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [7.4, 7.4]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")


jpl.plot_xy_data(nozzles=denotions, values_x=dnevi, values_y=Q_up_zrak, 
                 x_label=r"Dan", y_label=r"Toplotni tok uparjalnika [kW]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [6.4, 6.4]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")


jpl.plot_xy_data(nozzles=denotions, values_x=dnevi, values_y=P_komp_zrak, 
                 x_label=r"Dan", y_label=r"Moč kompresorja [kW]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [1, 1]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")


#COP_voda

#zanemarimo tlačne padce, izračunamo stanje vode
h_v_1 = PropsSI("H", "P", p_okolice, "T", T_podtalnice + 273, "Water") * 10**-3
h_v_2 = PropsSI("H", "P", p_okolice, "T", T_podtalnice - 3 + 273, "Water") * 10**-3
den_v = PropsSI("D", "P", p_okolice, "T", T_podtalnice - 3 + 273, "Water")

def moč_črpalke(Vf, density, g, H, f_d, mf, D, ni):
    a = density*g*H
    b = (8*mf**2)/(density*np.pi*D**5)
    c = H*f_d*b
    
    P_č = (Vf*(a + c))/ni
    return P_č

Q_up_voda = {}
Q_kond_voda = {}
P_komp_voda = {}
P_delo_voda = {}
COP_voda = {}
dnevi = {}
j = 0
for i in list(temp.keys()):
    Q_up = masni_pretok_voda[i] * (h_voda[i][0] - h_voda[i][-1])
    Q_kond = - masni_pretok_voda[i] * (h_voda[i][2] - h_voda[i][1])
    P_komp = masni_pretok_voda[i] * (h_voda[i][1] - h_voda[i][0])
    
    
    mf = - Q_up/(h_v_2 - h_v_1)
    vf = mf/den_v
    
    P_črp = moč_črpalke(vf, den_v, g, L_cevi, f_D, mf, D_cevi, izkoristek_črpalke) * 10**-3
    
    COP = Q_kond / (P_komp + P_črp)
    j+=1
    
    if i in ogrevani_dnevi:
        dnevi.update({i : j})
        P_komp_voda.update({i : P_komp})
        Q_up_voda.update({i : Q_up})
        Q_kond_voda.update({i : Q_kond})
        P_delo_voda.update({i : P_komp + P_črp})
        COP_voda.update({i : COP})
        
    else:
        n = f"negreje-{j}"
        dnevi.update({n : j})
        P_komp_voda.update({n : P_komp})
        Q_up_voda.update({n : Q_up})
        Q_kond_voda.update({n : Q_kond})
        P_delo_voda.update({n : P_komp + P_črp})

        COP_voda.update({n : COP})


denotions = ["day", "negreje"]
c_pal = {"day": "red", "negreje" : "blue"}
s_pal = {"day": ".", "negreje" : "X"}


jpl.plot_xy_data(nozzles=names, values_x=dnevi, values_y=COP_voda, 
                 x_label=r"Dan", y_label=r"COP [/]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [mejna_temp, mejna_temp]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")



jpl.plot_xy_data(nozzles=denotions, values_x=dnevi, values_y = Q_kond_voda, 
                 x_label=r"Dan", y_label=r"Toplotni tok kondenzatorja [kW]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [7.4, 7.4]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")


jpl.plot_xy_data(nozzles=denotions, values_x=dnevi, values_y=Q_up_voda, 
                 x_label=r"Dan", y_label=r"Toplotni tok uparjalnika [kW]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [6.4, 6.4]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")


jpl.plot_xy_data(nozzles=denotions, values_x=dnevi, values_y=P_delo_voda, 
                 x_label=r"Dan", y_label=r"Moč kompresor + črpalke [kW]", title="", 
                 figure_size = (13.8599/2.54, 11.0091/2.54),
                 scale="linear", marker_size = 50, 
                 show_grid = True, connect_line_color = {"day" : "red"},
                 point_opacity = 1, connect_points = None, connect_line_width = 1,
                 minor_ticks=True, color_palett = c_pal, symbol_palett = s_pal, 
                 legend_fontsize=9, legend_location="upper left", legend_param=(0,1),
                 x_tick_axis_spacing = 1, y_tick_axis_spacing = 1, plot_curve = False,
                 curve_points = [[[1, 365], [1, 1]]], curve_color = "black", curve_width = 2,
                 axis_fontsize=9, axis_scalar_format=False, show_legend=False, legend_marker_size = 100,
                 plot_font = "Palatino Linotype")


#Izračn sezonskih COP-jev

COP_zrak_delujoč = []
COP_voda_delujoč = []
for i in list(COP_zrak.keys()):
    if "negreje" not in i:
        COP_zrak_delujoč.append(COP_zrak[i])
        COP_voda_delujoč.append(COP_voda[i])


COP_sez_zrak = np.average(np.array(COP_zrak_delujoč))
COP_sez_voda = np.average(np.array(COP_voda_delujoč))

print("Sezonski COP za TČ zrak-voda znaša: " + f"{COP_sez_zrak}")
print("Sezonski COP za TČ voda-voda znaša: " + f"{COP_sez_voda}")








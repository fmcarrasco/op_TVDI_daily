# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 12:21:02 2017

@author: fcarrasco

Este script reune funciones para graficar analisis hechos
con los datos de MODIS diarios (EVI-LST-TVDI)
Productos en reticulas de 1km para

Temp. Superficial y EVI --> Calculo TVDI

"""
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from affine import Affine

import numpy as np
import pandas as pd
import datetime as dt
from scipy import stats

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sbn


def plot_triangle(evi, lst, tvdi, opc, zona):
    """
    Grafico de dispersion entre EVI y LST
    - Buscar que Parametros usar de acuerdo al dia Juliano
    - Graficar EVi vs LST
    """
    from modis_func import calc_day_tvdi
    from modis_func import param_tvdi
    # Calculate what parameter for TVDI need in terms of date of the file
    jul = calc_day_tvdi(opc)
    m, n, TSm = param_tvdi(jul, zona, opc)
    # ########
    x = evi.ReadAsArray().flatten()
    y = lst.ReadAsArray().flatten()
    z = tvdi.ReadAsArray().flatten()
    # Revisar este filtrado
    cond = np.logical_or(x == -3000, y == -3000)
    x[cond] = np.nan
    y[cond] = np.nan
    # Calculamos el TVDI
    # #################################################################
    # Graficamos la figura
    sbn.set(style='ticks', rc={'axes.grid': True,
                               'grid.linestyle':'--',
                               'grid.color':'#bac1ba'},
            palette='muted', color_codes=True, font_scale=0.6)

    fig, ax = plt.subplots(nrows=1, ncols=1, facecolor='white')
    #######################################
    ### Ploteamos en los distintos valores
    #######################################
    cnd = z < opc['itvdi'][0]
    ax.plot(x[cnd], y[cnd], color='#00734d', linestyle='None',\
            marker='.', ms=2.3, label='<0')
    cnd = np.logical_and(z >= opc['itvdi'][0], z < opc['itvdi'][1])
    ax.plot(x[cnd], y[cnd], color='#00734d', linestyle='None',\
            marker='.', ms=2.3, label='[0, 0.1)')
    cnd = np.logical_and(z >= opc['itvdi'][1], z < opc['itvdi'][2])
    ax.plot(x[cnd], y[cnd], color='#00e6a9', linestyle='None',\
            marker='.', ms=2.3, label='[0.1, 0.2)')
    cnd = np.logical_and(z >= opc['itvdi'][2], z < opc['itvdi'][3])
    ax.plot(x[cnd], y[cnd], color='#EADCDF', linestyle='None',\
            marker='.', ms=2.3, label='[0.2, 0.6)')
    cnd = np.logical_and(z >= opc['itvdi'][3], z < opc['itvdi'][4])
    ax.plot(x[cnd], y[cnd], color='#dfc27d', linestyle='None',\
            marker='.', ms=2.3, label='[0.6, 0.8)')
    cnd = z > opc['itvdi'][4]
    ax.plot(x[cnd], y[cnd], color='#dfc27d', linestyle='None',\
            marker='.', ms=2.3, label='> 0.8')
    # ###########################################################
    # Plot Limites
    x_lim = np.arange(-0.2, 1.1, 0.1)
    y_lim = m*x_lim + n
    ax.plot(x_lim, y_lim, color='#ce2727', linestyle='-.', alpha=0.7,\
            linewidth=1.5, label=u'Límite seco')
    ax.plot([-0.2, 1], [TSm, TSm], color='#3aaf42', alpha = 0.7,\
            linestyle='-.', label=u'Límite húmedo')
    # ###########################################################
    # Acomodamos Eje X
    ax.set_xlabel('EVI', fontsize=6)
    ax.set_xlim([-0.2, 1.1])
    # Acomodamos Eje y
    ax.set_ylabel('TS', fontsize=6)
    ax.set_ylim([280, 340])
    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator(2))
    ax.grid(b=True, which='minor', color='#bac1ba', linewidth=0.5)
    # Titulos
    yearjul = dt.datetime.strptime(opc['jul_f'], '%Y%j').strftime('%Y-%j')
    ax.set_title('Zona: ' + zona + u'; Año: ' + yearjul, fontsize=11)
    # Colocamos leyenda
    lgnd = ax.legend(loc='best', fancybox=True, prop={'size': 11}, handletextpad=0.1)
    # Cambiamos el tamagno de marcadores en leyenda
    lgnd.legendHandles[0]._legmarker.set_markersize(10)
    lgnd.legendHandles[1]._legmarker.set_markersize(10)
    lgnd.legendHandles[2]._legmarker.set_markersize(10)
    lgnd.legendHandles[3]._legmarker.set_markersize(10)
    lgnd.legendHandles[4]._legmarker.set_markersize(10)
    lgnd.legendHandles[5]._legmarker.set_markersize(10)
    lgnd.legendHandles[6]._legmarker.set_markersize(10)
    lgnd.get_frame().set_facecolor('white')
    # Guardamos la figura
    sbn.despine(ax=ax, offset=10, trim=True)
    fechahoy = dt.datetime.today().strftime('%Y%j%H%M')
    figname = opc['sf'] + zona + '_' + yearjul +\
              '_' + fechahoy + '_diagnostico.jpg'
    plt.savefig(figname, dpi=200)
    plt.close(fig)

def CalcTotalPoints(region, opc):
    """
    """
    from modis_func import get_lim_mask

    zcfolder = opc['zc']
    # ########################
    if 'NEA' in region:
        mascara = gdal.Open(zcfolder+'MASK_NEA')
    elif 'NORTE' in region:
        mascara = gdal.Open(zcfolder+'MASK_NORTE')
    elif 'PS' in region:
        mascara = gdal.Open(zcfolder+'MASK_PS')
    elif 'PN' in region:
        mascara = gdal.Open(zcfolder+'MASK_PN')
    elif 'NN' in region:
        mascara = gdal.Open(zcfolder+'MASK_NN')
    else:
        import sys
        print('------------------------------------------')
        print('Region desconocida, el valor: ' + region)
        print('No corresponde a ninguna de las utilizadas')
        print('actualmente. Revisar parametro ZONA')
        print('------------------------------------------')
        sys.exit()
    # #####################################################
    mk = mascara.ReadAsArray()
    x1, x2, y1, y2 = get_lim_mask(mk)
    Xoff = int(x1)
    Yoff = int(y1)
    Cols = int(x2 - x1)
    Rows = int(y2 - y1)
    mask_part = mascara.ReadAsArray(Xoff, Yoff, Cols, Rows)
    n_real_data = np.sum(mask_part == 1)

    return n_real_data

def ClassTVDI(f_tvdi, opc, clas):
    """
    Genera un array con la clasificacion de los datos segun el
    array clas
    """
    dato = f_tvdi.ReadAsArray()
    dato[dato == -3000] = np.nan
    out_arr = []
    a1 = np.nansum(dato < clas[0])  # Primer dato de clas
    out_arr.append(a1)
    for low, up in zip(clas, clas[1:]):
        a = np.nansum( np.logical_and(dato >= low, dato < up) )
        out_arr.append(a)
    a2 = np.nansum(dato >= clas[-1])  # Ultimo dato de clas
    out_arr.append(a2)

    return out_arr


def CalcTablaResumen(NN, PN, PS, opc):
    """
    Generamos una tabla con el resumen sobre los datos
    """
    # Calculamos el total de puntos en cada zona
    pts_NN = CalcTotalPoints('NN', opc)
    pts_PN = CalcTotalPoints('PN', opc)
    pts_PS = CalcTotalPoints('PS', opc)
    # Clasificamos
    arr_clas = [0., 0.1, 0.2, 0.6, 0.8, 1.]
    clas_NN = ClassTVDI(NN, opc, arr_clas)
    clas_PN = ClassTVDI(PN, opc, arr_clas)
    clas_PS = ClassTVDI(PS, opc, arr_clas)
    # Calculamos porcentajes
    porc_NN = [100.*(item/pts_NN) for item in clas_NN]
    porc_PN = [100.*(item/pts_PN) for item in clas_PN]
    porc_PS = [100.*(item/pts_PS) for item in clas_PS]
    #
    falt_NN = 100*(pts_NN - sum(clas_NN))/pts_NN
    falt_PN = 100*(pts_PN - sum(clas_PN))/pts_PN
    falt_PS = 100*(pts_PS - sum(clas_PS))/pts_PS
    #
    porc_NN.append(falt_NN)
    porc_PN.append(falt_PN)
    porc_PS.append(falt_PS)
    # Generamos un Pandas Dataframe para luego guardar como TXT y EXCEL
    a_col = ['TVDI < 0.0', 'TVDI en [0.0, 0.1)', 'TVDI en [0.1, 0.2)',\
             'TVDI en [0.2, 0.6)', 'TVDI en [0.6, 0.8)', 'TVDI en [0.8, 1.0)',\
             'TVDI >= 1', 'Faltantes']
    df = pd.DataFrame(index=np.arange(0,8))
    df = df.assign(Clas=a_col)
    df = df.assign(NN=porc_NN)
    df = df.assign(PN=porc_PN)
    df = df.assign(PS=porc_PS)
    #
    df = df.apply(pd.to_numeric, errors='ignore')
    # Guardamos el archivo
    fechahoy = dt.datetime.today().strftime('%Y%j%H%M')
    yearjul = dt.datetime.strptime(opc['jul_f'], '%Y%j').strftime('%Y%j')
    nombre = 'Resumen_TVDI_' + yearjul + '_' + fechahoy
    df.to_excel(opc['sf'] + nombre + '.xlsx', float_format='%.4f',
                index=False, na_rep='9999', encoding='ISO-8859-1')
    df.to_csv(opc['sf'] + nombre + '.csv', sep=';', decimal=',',
              float_format='%.4f', index=False)

    return clas_NN

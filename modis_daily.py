# -*- coding: utf-8 -*-
"""
Created on Wed 20/02 2017

@author: fcarrasco

Este script reune funciones para trabajar con los
datos HDF de MODIS.
Esta adaptado para obtener un TIFF en resolucion de 1km
de TVDI a partir de datos diarios.

Hay distintos pasos para esta tarea:
1) Estimar EVI a partir de datos de bandas de reflectividad
2) Leer datos de TS
3) Combinar ambos datos previos y utilizar las mismas zonas
  y parametros definidos para estimar TVDI en periodos de 16 dias.
4) Chequear que los valores del espacio EVI-TS esten dentro del rango
de los triangulos de cada zona.

Los codigos escritos fueron inspirados en los
codigos publicados en:

https://jgomezdans.github.io/gdal_notes/reprojection.html

"""
# Importing packages

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from affine import Affine
import numpy as np
import pandas as pd
import datetime as dt

from modis_func import extract_fname
from modis_func import getVar_from_HDF
from modis_func import transform
from modis_func import filtrar
from modis_func import get_bit_error
from modis_func import reproyectar
from modis_func import subset
from modis_func import get_lim_mask
from modis_func import modis_to_latlon
from modis_func import CalcEVI


def WrkWithLST_files(l_files, opci):
    """
    Main function to work with  4 LST files (each one in 1km resolution)
    -Extract data from filename
    -Extract data from each HDF to a GTIFF with Float format
    -Filter data using QC variable in HDF
    -Make a mosaic with the four LST filtered images (Save a GTIFF)
    -Extract LST data to each of the three zones that ORA monitor
    -Return the three complete path for each file
    """
    # PASO 1
    file_data = extract_fname(l_files[0])
    opci['f_data'] = file_data
    # PASO 2
    print('Procesando datos de LST para: ' +\
          file_data['fecha'].strftime('%Y-%j'))
    opci['l_file'].write('Procesando datos de EVI para: ' +\
                         file_data['fecha'].strftime('%Y-%j') + '\n')
    # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    var = getVar_from_HDF(l_files, 'LST_Day_1km')
    qc = getVar_from_HDF(l_files, 'QC_Day')
    data = []
    for idx, item in enumerate(var):
        # Que archivo trabajamos
        opci['l_file'].write(l_files[idx] + '\n')
        # Extract data from file
        ds = gdal.Open(item)
        NX = ds.RasterXSize
        NY = ds.RasterYSize
        fill_value = np.float64(ds.GetMetadata()['_FillValue'])
        scale_factor = np.float64(ds.GetMetadata()['scale_factor'])
        # Transform; Filter and Composite data
        data = transform(ds.ReadAsArray(), fill_value, scale_factor,
                         opci['l_file'])
        ds_f = gdal.Open(qc[idx])
        data_flt = filtrar(data, ds_f.ReadAsArray(), opci)
        if idx == 0:
            auxdata = np.zeros((NX, NY,4))
            auxdata[:] = np.NAN
        auxdata[:,:,idx] = data_flt
    m_composite = np.empty((2*NX, 2*NY))
    m_composite[0:NX,0:NY] = auxdata[:,:,0]
    m_composite[NX:,0:NY] = auxdata[:,:,1]
    m_composite[0:NX,NY:] = auxdata[:,:,2]
    m_composite[NX:,NY:] = auxdata[:,:,3]
    # Liberamos memoria
    del auxdata, data, data_flt
    # Reproyect
    mos_lst = reproyectar(var[0], m_composite, 'sinu', 'LST_Day_', opci)
    # Subset in each of the zones
    mos_lst_sub0 = subset(mos_lst, 'LST_Day', 'NN', opci)
    mos_lst_sub1 = subset(mos_lst, 'LST_Day', 'PN', opci)
    mos_lst_sub2 = subset(mos_lst, 'LST_Day', 'PS', opci)

    return mos_lst_sub0, mos_lst_sub1, mos_lst_sub2


def WrkWithREFLECT_files(l_files, refl, opc):
    """
    Main function to work with  4 LST files (each one in 1km resolution)
    -Extract data from filename
    -Extract data from each HDF to a GTIFF with Float format
    -Filter data using QC variable in HDF
    -Make a mosaic with the four LST filtered images (Save a GTIFF)
    -Extract LST data to each of the three zones that ORA monitor
    -Return the three complete path for each file
    """
    # PASO 1
    file_data = extract_fname(l_files[0])
    opc['f_data'] = file_data
    # PASO 2
    print('Procesando datos de Reflectividad ' + refl + ' para: ' +\
          file_data['fecha'].strftime('%Y-%j'))
    opc['l_file'].write('Procesando datos de EVI para: ' +\
                        file_data['fecha'].strftime('%Y-%j') + '\n')
    # # # # # # # # # # # # # # # # # # # #
    prf = ['QC_500m_1', 'sur_refl_b01_1', 'sur_refl_b02_1', 'sur_refl_b03_1']
    qc = getVar_from_HDF(l_files, 'QC_500m_1')
    var = getVar_from_HDF(l_files, refl)
    # # # # # # # # #
    for idx, item in enumerate(var):
        # Que archivo trabajamos
        opc['l_file'].write(l_files[idx] + '\n')
        # Extract data from file
        ds = gdal.Open(item)
        NX = ds.RasterXSize
        NY = ds.RasterYSize
        fill_value = np.float64(ds.GetMetadata()['_FillValue'])
        scale_factor = np.float64(0.0001)
        # Transform; Filter and Composite data
        data = transform(ds.ReadAsArray(), fill_value, scale_factor,
                         opc['l_file'])
        ds_f = gdal.Open(qc[idx])
        valores = np.unique(ds_f.ReadAsArray())
        if idx == 0:
            auxdata = np.zeros((NX, NY,4))
            auxdata[:] = np.NAN
        auxdata[:,:,idx] = data
    # ### END of LOOP ###
    m_composite = np.empty((2*NX, 2*NY))
    m_composite[0:NX,0:NY] = auxdata[:,:,0]
    m_composite[NX:,0:NY] = auxdata[:,:,1]
    m_composite[0:NX,NY:] = auxdata[:,:,2]
    m_composite[NX:,NY:] = auxdata[:,:,3]
    # Liberamos memoria
    del auxdata, data
    #print(refl_item[0])
    mos_lst = reproyectar(var[0], m_composite, 'sinu',
                          refl + '_', opc)
    # Subset in each of the zones
    mos_lst_sub0 = subset(mos_lst, refl, 'NN', opc)
    mos_lst_sub1 = subset(mos_lst, refl, 'PN', opc)
    mos_lst_sub2 = subset(mos_lst, refl, 'PS', opc)


    return mos_lst_sub0, mos_lst_sub1, mos_lst_sub2


if __name__ == "__main__":

    import glob
    import subprocess
    import time
    import os
    import matplotlib
    import matplotlib.pyplot as plt
    from modis_plt import plot_triangle
    from modis_func import calc_day_tvdi
    from modis_func import CalcTVDI
    from modis_plt import CalcTablaResumen
    # ######################################################################
    # Parametros HDF
    path = 'e:/TVDI/hdf/daily/'
    year = '2019'
    month = '01'
    day = '20'
    c_date = year + '-' + month + '-' + day
    vfile = ['MYD09GA', 'MYD11A1']
    zona = 'NEA'  # Subset: 'NN'; 'NEA'; 'NORTE';'PN';'PS'
    p_path = './parametros_csv/'  # Carpeta Parametros
    z_path = 'c:/Felix/ORA/MODIS/Zonas_climaticas/'  # Carpeta Zonas Clim HDR
    fechahoy = dt.datetime.today().strftime('%Y%j')
    juld     = dt.datetime.strptime(c_date, '%Y-%m-%d').strftime('%Y%j')
    outfolder = './salida_modis_diaria/'
    tmpfolder = outfolder + 'tmp/'
    salidas = outfolder + 'tvdi/tvdi_' + year + month + day +\
              '_' + fechahoy + '/'
    # Opciones Filtrado
    opt_filtrar = 1  # 0: No Filtrar; 1: MB y B; 2: MB
    # Opciones en salidas y Reproyecciones:
    px_sz = 1000  # Size of pixel in meters
    out_format = 'GTiff'
    # Valores a graficar
    i_tvdi = [0., 0.1, 0.2, 0.6, 0.8, 1.]
    # #####################################################################
    # Archivo log
    logfile = 'tvdi_diario_' + juld + '_' + fechahoy + '.log'
    f = open(outfolder + logfile, 'w')
    f.write('Calculo de TVDI con datos DIARIOS\n')
    f.write('\n')
    f.write('Datos HDF en: '+ path + '/\n')
    f.write('Datos zonas climaticas en: ' + z_path + '\n')
    f.write('Salidas temporales en: ' + tmpfolder + '\n')
    f.write('Salidas TVDI en: ' + salidas + '\n')
    # #####################################################################
    # Start working with daily temperature
    os.makedirs(salidas, exist_ok=True)
    os.makedirs(tmpfolder, exist_ok=True)
    y_path = path + year + month + day + '/'
    prefix = vfile[1] + '*.hdf'
    file_tmp = glob.glob(y_path + prefix)
    opc = {'p_csv':p_path, 'zc':z_path, 'of':outfolder, 'tmpf': tmpfolder,
           'sf':salidas, 'filt':opt_filtrar, 'jul_f':juld,
           'l_file': f,
           'fmt_out':out_format, 'px_sz':px_sz,
           'itvdi':i_tvdi}

    lst_NN, lst_PN, lst_PS = WrkWithLST_files(file_tmp, opc)
    # #####################################################################
    # Start Working with data Reflectivity
    y_path = path + year + month + day + '/'
    prefix = vfile[0] + '*.hdf'
    file_ref = glob.glob(y_path + prefix)
    #'sur_refl_b01_1', 'sur_refl_b02_1', 'sur_refl_b03_1'
    rf1_NN, rf1_PN, rf1_PS = WrkWithREFLECT_files(file_ref, 'sur_refl_b01_1',
                                                  opc)
    rf2_NN, rf2_PN, rf2_PS = WrkWithREFLECT_files(file_ref, 'sur_refl_b02_1',
                                                  opc)
    rf3_NN, rf3_PN, rf3_PS = WrkWithREFLECT_files(file_ref, 'sur_refl_b03_1',
                                                  opc)

    # Calculate EVI from Reflectivity
    # #####################################################################
    print('EVI zona NN')
    evi_NN = CalcEVI(rf1_NN, rf2_NN, rf3_NN, opc)
    print('EVI zona PN')
    evi_PN = CalcEVI(rf1_PN, rf2_PN, rf3_PN, opc)
    print('EVI zona PS')
    evi_PS = CalcEVI(rf1_PS, rf2_PS, rf3_PS, opc)
    # #####################################################################
    # Calculamos TVDI
    tvdi_NN = CalcTVDI(evi_NN, lst_NN, opc, 'NN')
    tvdi_PN = CalcTVDI(evi_PN, lst_PN, opc, 'PN')
    tvdi_PS = CalcTVDI(evi_PS, lst_PS, opc, 'PS')
    # Graficamos plano EVI vs LST
    plot_triangle(evi_NN, lst_NN, opc, 'NN')
    plot_triangle(evi_PN, lst_PN, opc, 'PN')
    plot_triangle(evi_PS, lst_PS, opc, 'PS')
    # Generamos tabla con resumen de datos
    CalcTablaResumen(tvdi_NN, tvdi_PN, tvdi_PS, opc)
    # Borramos memoria
    lst_NN = None
    lst_PN = None
    lst_PS = None
    rf1_NN = None
    rf1_PN = None
    rf1_PS = None
    rf2_NN = None
    rf2_PN = None
    rf2_PS = None
    rf3_NN = None
    rf3_PN = None
    rf3_PS = None
    f.close()

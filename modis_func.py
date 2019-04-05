"""
Function to use for calculation of daily TVDI

"""

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from affine import Affine
import numpy as np
import pandas as pd
import datetime as dt

def extract_fname(filename):
    """
    Funcion para extraer las distintas partes del nombre
    de archivo HDF para MODIS segun la convencion descrita
    en el siguiente archivo:

    https://docs.google.com/document/d/
    1qSTBaasz3-4wOb9VSOXm3rDLKKOgL1b_lcge-4mOkyk/edit?usp=sharing

    """
    aux = filename.split('.')
    short_name = aux[0].split('\\')[1]
    data_date =  dt.datetime.strptime(aux[1], 'A%Y%j')
    tile = aux[2]
    version = aux[3]
    create_date =  dt.datetime.strptime(aux[4], '%Y%j%H%M%S')
    dictionary = {'nvar': short_name, 'fecha': data_date,
                  'tile': tile, 'version': version,
                  'fecha_c': create_date}
    return dictionary


def getVar_from_HDF(archivos, prefix):
    """
    Funcion para extraer los datasets necesarios
    de un listado de archivos HDF.
    Conviene revisar el archivo previamente para
    colocar bien el prefijo
    """
    listado = []
    for arc in archivos:
        ds = gdal.Open(arc)
        for sd, descr in ds.GetSubDatasets():
            if (prefix in descr.split(' ')):
                listado.append(sd)
        ds = None

    return listado


def transform(M, fv, sf, f):
    """
    Transforma datos del HDF a valores punto flotante y utiliza la
    transformacion:

    dato = dato_orig * scale_factor
    where dato == _FillValue --> NaN
    M: dato_orig
    fv: _FillValue
    sf: scale_factor
    """
    arr_out = np.zeros(np.shape(M))
    arr_out = np.float64(M)*sf
    arr_out[M == fv] = np.nan
    p_total = np.shape(M)[0]*np.shape(M)[1]
    filtered = np.sum(M == fv)
    porc = (np.float64(filtered)/p_total)*100.
    porc2 = "%.4f" % porc
    f.write('Puntos Totales = ' +str(p_total) +\
            ' Puntos filtrados = ' + str(filtered) +\
            ' [porcentual: ' + porc2 + ']\n')

    return arr_out


def filtrar(M, flt, opc):
    """
    Filtra datos a NaN de la matriz M, usando los datos de la matriz flt

    Notar que los datos en flt, hay que pasarlos a bit y luego usarlos.

    opt = 1 -> Deja datos buenos y muy buenos
    opt = 2 -> Deja datos muy buenos.
    """
    values = np.unique(flt)
    #print(values)
    bad_values = get_bit_error(values, opc['filt'])
    #print(bad_values)
    M1 = M
    t_filt = 0.
    for item in bad_values:
        M1[flt == item] = np.nan
        t_filt = t_filt + np.sum(flt == item)
    p_total = np.shape(M)[0]*np.shape(M)[1]
    porc = (np.float64(t_filt)/p_total)*100.
    porc2 = "%.4f" % porc
    opc['l_file'].write('Puntos Totales = ' +str(p_total) +\
                        ' Puntos filtrados = ' + str(t_filt) +\
                        ' [porcentual: ' + porc2 + ']\n')

    return M1


def get_bit_error(array, opt):
    """
    Funcion que para cada item en array, lo convierte a
    formato binario para poder calcular cuales son los
    que permiten reconocer errores
    """
    bad_values = []
    for item in array:
        bin_val = str(bin(int(item)))[2:].zfill(8)
        bit_01 = str(bin_val)[0:2]
        if opt == 1:
            if ('10' in bit_01) or ('11' in bit_01):
                bad_values.append(item)
        elif opt == 2:
            if ('10' in bit_01) or ('11' in bit_01) or ('01' in bit_01):
                bad_values.append(item)

    return bad_values


def reproyectar(dataset, M, out_proj, prefix, opc):
    """
    Esta funcion reproyecta un dataset a la coordenada espacial epsg_to
    y con el pixel deseado (pixel_size).
    Se guardan en el formato especificado.

    """
    g = gdal.Open ( dataset )
    # Set up the two Spatial Reference systems.
    in_proj = osr.SpatialReference(wkt=g.GetProjection())
    # Open the original dataset, and get the geotransform
    # Get the Geotransform vector
    geo_t = g.GetGeoTransform ()
    InCols = np.shape(M)[0]  # Raster xsize
    InRows = np.shape(M)[1]  # Raster ysize
    # Name for New file_shp
    fhoy = dt.datetime.today().strftime('%Y%j')
    in_file_name = opc['tmpf'] + prefix + '_' + fhoy + '_orig.tiff'
    out_file_name = opc['tmpf'] + prefix + '_' + fhoy + '_repr.tiff'
    # Create an in-memory raster dataset
    mem_drv = gdal.GetDriverByName( opc['fmt_out'] )
    InDS = mem_drv.Create(in_file_name, InCols, InRows, 1, gdal.GDT_Float64 )
    InDS.SetGeoTransform( geo_t )
    InDS.SetProjection ( in_proj.ExportToWkt() )
    InDS.GetRasterBand(1).WriteArray(M)
    InDS.FlushCache()
    # Reproyectamos la imagen al tamagno de pixel deseado:
    if out_proj == 'sinu':
        gdal.Warp(out_file_name, in_file_name,
                  dstSRS=in_proj.ExportToWkt(),
                  dstNodata=-3000.,
                  xRes=1000, yRes=1000)
    elif out_proj == 'latlon':
        gdal.Warp(out_file_name, in_file_name,
                  dstSRS='EPSG:4326',
                  dstNodata=-3000.)
    # Lineas para el LOGFILE
    opc['l_file'].write('\n')
    opc['l_file'].write('Archivo tiff original: ' +\
                        in_file_name +'\n')
    opc['l_file'].write('Archivo tiff reproyectado: ' +\
                         out_file_name +'\n')
    opc['l_file'].write('\n')
    # Devolvemos el archivo de salida
    return out_file_name


def subset(dataset, prefix, region, opc):
    """
    Realizar un subset sobre el dataset de input.
    Hay 4 regiones establecidas:
    -------------------------
    prefijo --> archivo
    NEA --> MASK_NEA
    NORTE --> MASK_NORTE
    PN --> MASK_PN
    PS --> MASK_PS
    NN --> MASK_NN
    --------------------------
    opc: Diccionario con las opciones donde guardar y encontrar
    los respectivos archivos.
    """
    folder = opc['tmpf']
    zcfolder = opc['zc']
    fecha = opc['f_data']['fecha'].strftime('%Y%j')
    if 'NEA' in region:
        mk_file = zcfolder+'MASK_NEA'
        mascara = gdal.Open(mk_file)
    elif 'NORTE' in region:
        mk_file = zcfolder+'MASK_NORTE'
        mascara = gdal.Open(mk_file)
    elif 'PS' in region:
        mk_file = zcfolder+'MASK_PS'
        mascara = gdal.Open(mk_file)
    elif 'PN' in region:
        mk_file = zcfolder+'MASK_PN'
        mascara = gdal.Open(mk_file)
    elif 'NN' in region:
        mk_file = zcfolder+'MASK_NN'
        mascara = gdal.Open(mk_file)
    else:
        import sys
        print('------------------------------------------')
        print('Region desconocida, el valor: ' + region)
        print('No corresponde a ninguna de las utilizadas')
        print('actualmente. Revisar parametro region')
        print('------------------------------------------')
        sys.exit()
    data_2 = modis_to_latlon(gdal.Open(dataset), mascara, prefix, region, opc)
    a  = Affine.from_gdal(*mascara.GetGeoTransform())
    mk = mascara.ReadAsArray()
    ####################################
    # Obtenemos los limites de la mascara
    # y leemos en esos limites
    x1, x2, y1, y2 = get_lim_mask(mk)
    Xoff = int(x1)
    Yoff = int(y1)
    Cols = int(x2 - x1)
    Rows = int(y2 - y1)
    data_part = data_2.ReadAsArray(Xoff, Yoff, Cols, Rows)
    mask_part = mascara.ReadAsArray(Xoff, Yoff, Cols, Rows)
    ###################################################
    # Lon-Lat del punto de subset en la esquina Sup-Izq
    lon1, lat1 = a*(x1, y1)
    # Enmascaramos los datos
    data_part[mask_part == 0] = np.NaN
    n_real_data = np.sum(mask_part == 1)
    # Abrimos un archivo para guardar
    mem_drv = gdal.GetDriverByName( 'Gtiff' )
    fnames = folder + prefix + '_' + region + '_' + fecha + '_subset.tiff'
    dest = mem_drv.Create(fnames, Cols, Rows, 1, gdal.GDT_Float64)
    opc['l_file'].write('\n')
    opc['l_file'].write('Subset en: ' + fnames + '\n')
    opc['l_file'].write('\n')
    # Datos Para Guardar
    gt = mascara.GetGeoTransform()
    new_geo = [lon1, gt[1], gt[2], lat1, gt[4], gt[5]]
    dest.SetGeoTransform( new_geo )
    dest.SetProjection ( mascara.GetProjection() )
    data_part[np.isnan(data_part)] = -3000.
    dest.GetRasterBand(1).SetNoDataValue(-3000)
    dest.GetRasterBand(1).WriteArray(data_part)
    dest.FlushCache()

    return dest


def get_lim_mask(mascara):
    """
    Funcion para determinar los limites superior e inferior
    de la mascara para poder definir los limites del
    recuadro enmascarado.

    Input:
    mascara: matriz con 0 y 1 donde 1 representa los valores
    a considerar.

    Output:
    x1, x2, y1, y2 : Puntos a definir

    """
    ########################################################
    a1 = []
    a1 = np.argmax(mascara, axis=0)
    y1 = min(a1[np.nonzero(a1)]) - 3
    ########################################################
    a1 = []
    a1 = np.argmax(np.flipud(mascara), axis=0)
    y2 = mascara.shape[0] - min(a1[np.nonzero(a1)]) - 1 + 3
    ########################################################
    a1 = []
    a1 = np.argmax(mascara, axis=1)
    x1 = min(a1[np.nonzero(a1)]) - 3
    ########################################################
    a1 = []
    a1 = np.argmax(np.fliplr(mascara), axis=1)
    x2 = mascara.shape[1] - min(a1[np.nonzero(a1)]) - 1 + 3
    ########################################################
    return x1, x2, y1, y2


def modis_to_latlon(g, mask, region, prefix, opc):
    """
    Funcion para pasar de la proyeccion de MODIS a la
    proyeccion Lat-Lon con Datum WGS84.
    """
    ##################################################
    # Definimos las proyecciones
    # y definimos la funcion para transformar entre coordenadas
    ##################################################
    wkt_from = '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'
    epsg_to = 4326
    wgs84 = osr.SpatialReference ()
    wgs84.ImportFromEPSG ( epsg_to )
    modis = osr.SpatialReference ()
    modis.ImportFromProj4( wkt_from )
    tx = osr.CoordinateTransformation ( modis, wgs84 )
    ##################################################
    # Abrimos el dato MODIS que viene ya como GdalDataSet
    # Vector de coordenadas Geograficas y tamagno imagen
    ##################################################
    geo_t = g.GetGeoTransform ()
    x_size = g.RasterXSize # Raster xsize
    y_size = g.RasterYSize # Raster ysize
    #################################################
    # Trabajamos los bordes de la imagen proyectada en WGS84
    #################################################
    (ulx, uly, ulz ) = tx.TransformPoint( geo_t[0], geo_t[3])
    (lrx, lry, lrz ) = tx.TransformPoint( geo_t[0] + geo_t[1]*x_size, \
                                          geo_t[3] + geo_t[5]*y_size )
    #################################################
    # Archivo de mascara para extraer tamaño de pixel y bordes
    #################################################
    geotransform = mask.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    x_size_m = mask.RasterXSize
    y_size_m = mask.RasterYSize
    ##################################################
    # Abrimos el archivo temporal
    ##################################################
    mem_drv = gdal.GetDriverByName( 'GTiff' )
    x_px_space = pixelWidth
    y_px_space = -pixelHeight
    x_size_est = int((lrx - ulx)/x_px_space)
    y_size_est = int((uly - lry)/y_px_space)
    ##################################################
    # Creamos el archivo que guarda la
    ##################################################
    fecha = opc['f_data']['fecha'].strftime('%Y-%j')
    dest = mem_drv.Create(opc['tmpf'] + prefix + '_' +\
                          '_' + fecha + '.tiff',
                          x_size_m, y_size_m, 1, gdal.GDT_Float64)
    # Calcular las nuevas coordenadas geograficas
    new_geo = ( ulx, x_px_space, geo_t[2], \
                uly, geo_t[4], -y_px_space )
    # Set the geotransform
    dest.SetGeoTransform( geotransform )
    dest.SetProjection ( wgs84.ExportToWkt() )
    # Perform the projection/resampling
    gdal.ReprojectImage( g, dest, \
                         modis.ExportToWkt(), wgs84.ExportToWkt(), \
                         gdal.GRA_NearestNeighbour )
    dest.GetRasterBand(1).SetNoDataValue(-3000)

    return dest


def CalcEVI(rf1, rf2, rf3, opc):
    """
    Calculate EVI from three reflectivity bands.
    The formula is taken from
    Huete et al., 2002
    "Overview of the radiometric and biophysical
    performance of the MODIS vegetation indices"
                        NIR - RED
    EVI = G x ______________________________
               (NIR + C1*RED - C2*Blue + L)
    rf1: RED Band
    rf2: NIR Band
    rf3: Blue Band
    """
    # Define Parameters
    G = 2.5
    L = 1.
    C1 = 6.
    C2 = 7.5
    # Read bands as arrays
    red = rf1.ReadAsArray()
    red[red == -3000] = np.nan
    nir = rf2.ReadAsArray()
    nir[nir == -3000] = np.nan
    blue = rf3.ReadAsArray()
    blue[blue == -3000] = np.nan
    # Calculate numerador y denominador
    nume = G * np.where(np.isnan(nir),
                        np.nan, nir - red)
    nume[np.isnan(red)] = np.nan
    deno = nir + C1 * red - C2 * blue + L
    # Calculate EVI
    evi = np.empty(np.shape(red))
    evi[:] = np.nan
    evi = np.where(deno == 0, np.nan, np.divide(nume, deno))
    evi[np.isnan(evi)] = -3000
    # Save data as GeoTiff
    fecha = opc['f_data']['fecha'].strftime('%Y-%j')
    aux = rf1.GetDescription().split('_')  # Get file name and split by '_'
    suffix = '_'.join(aux[6::])  # Join by '_' using only the date and zone
    namefile = opc['tmpf'] + 'EVI_' + suffix
    mem_drv = gdal.GetDriverByName( 'GTiff' )
    dest = mem_drv.Create(namefile, rf1.RasterXSize, rf1.RasterYSize, 1,
                          gdal.GDT_Float64)
    dest.SetGeoTransform( rf1.GetGeoTransform() )
    dest.SetProjection ( rf1.GetProjection() )
    dest.GetRasterBand(1).SetNoDataValue(-3000)
    dest.GetRasterBand(1).WriteArray(evi)
    dest.FlushCache()

    return dest


def calc_day_tvdi(opc):
    """
    Calculamos que dia le corresponde entre los periodos de 16 dias
    que actualmente se hace seguimiento en la ORA.
    Con esto se define que parametros utilizar para calcular
    el TVDI.
    """
    datefile = dt.datetime.strptime(opc['jul_f'], '%Y%j')  # Fecha del archivo
    jul_file = datefile.timetuple().tm_yday
    # Dias que se hacen seguimiento.
    jevi = [9, 25, 41, 57, 73, 89, 281, 297, 313, 329, 345, 361]
    #
    if np.logical_and(jul_file >= jevi[0], jul_file < jevi[1]):
        jul_out = jevi[0]
    elif np.logical_and(jul_file >= jevi[1], jul_file < jevi[2]):
        jul_out = jevi[1]
    elif np.logical_and(jul_file >= jevi[2], jul_file < jevi[3]):
        jul_out = jevi[2]
    elif np.logical_and(jul_file >= jevi[3], jul_file < jevi[4]):
        jul_out = jevi[3]
    elif np.logical_and(jul_file >= jevi[4], jul_file < jevi[5]):
        jul_out = jevi[4]
    elif np.logical_and(jul_file >= jevi[5], jul_file < jevi[5] + 16):
        jul_out = jevi[5]
    elif np.logical_and(jul_file >= jevi[6], jul_file < jevi[7]):
        jul_out = jevi[6]
    elif np.logical_and(jul_file >= jevi[7], jul_file < jevi[8]):
        jul_out = jevi[7]
    elif np.logical_and(jul_file >= jevi[8], jul_file < jevi[9]):
        jul_out = jevi[8]
    elif np.logical_and(jul_file >= jevi[9], jul_file < jevi[10]):
        jul_out = jevi[9]
    elif np.logical_and(jul_file >= jevi[10], jul_file < jevi[11]):
        jul_out = jevi[10]
    elif np.logical_or(jul_file >= jevi[11], jul_file < jevi[1]):
        jul_out = jevi[11]
    else:
        print('Fecha sin seguimiento de TVDI --> NaN')
        jul_out = np.nan

    return jul_out


def param_tvdi(djul, zona, opc):
    """
    Devuelve los parametros del TVDI a partir del dia Juliano
    indicado y la zona.
    """

    filename = 'summary_file_' + zona + '.csv'
    df = pd.read_csv(opc['p_csv'] + filename, sep=';', decimal=',',
                     encoding='ISO-8859-1', index_col=0)

    ### Parametros TVDI ####
    m = df['m'][djul]
    n = df['n'][djul]
    TSm = df['TSmin'][djul]
    del df

    return m, n, TSm


def save_Gtiff(dato, t_file, zona, opc):
    """
    Funcion para guardar un GeoTiff en misma proyeccion:
    dato: Matriz Numpy con los datos a guardar
    t_file: Archivo Geotiff con las caracteristicas para guardar el dato
    DEBEN coincidir las dimensiones, si no, va a tirar error
    opc: listado de opciones con las carpetas donde guardar

    """
    # Abrimos un archivo para guardar
    mem_drv = gdal.GetDriverByName( 'Gtiff' )
    # Parametros Salida
    yearjul = dt.datetime.strptime(opc['jul_f'], '%Y%j').strftime('%Y%j')
    fechahoy = dt.datetime.today().strftime('%Y%j%H%M')
    outname = 'TVDI_' + yearjul + '_' + fechahoy + '_' + zona + '.tiff'
    outfolder = opc['sf']
    # Datos GeoTiff
    geo_t = t_file.GetGeoTransform ()
    x_size = t_file.RasterXSize # Raster xsize
    y_size = t_file.RasterYSize # Raster ysize
    opc['l_file'].write('Guardando archivo TVDI en: ' + outfolder +\
                        outname + '\n')
    opc['l_file'].write('\n')
    dest = mem_drv.Create(outfolder + outname, x_size, y_size, 1,
                          gdal.GDT_Float64)
    dest.SetGeoTransform( geo_t )
    dest.SetProjection ( t_file.GetProjection() )
    dato[np.isnan(dato)] = -3000.
    dest.GetRasterBand(1).SetNoDataValue(-3000)
    dest.GetRasterBand(1).WriteArray(dato)
    dest.FlushCache()

    return dest


def CalcTVDI(evi_f, lst_f, opc, zona):
    """
    Funcion para calcular TVDI.
    Se genera un GeoTIFF para la zona especificada
    a partir de los datos de EVI y LST.

    La formula del TVDI es la siguiente:
                  TS - TSmin
    TVDI = ------------------------
              m*EVI + n - TSmin

    Los parametros se calculan a partir de la fecha de los
    archivos
    """
    # Primero calculamos el dia Juliano correspondiente a los parametros.
    jul = calc_day_tvdi(opc)
    m, n, TSm = param_tvdi(jul, zona, opc)
    # Leemos las matrices de datos
    lst = lst_f.ReadAsArray()
    evi = evi_f.ReadAsArray()
    cond = np.logical_or(lst == -3000, evi == -3000)
    lst[cond] = np.nan
    evi[cond] = np.nan
    #
    nume = lst - TSm
    deno = m*evi + n - TSm
    nume[deno==0] = np.nan
    deno[deno==0] = np.nan
    tvdi = np.divide(nume, deno)
    tvdi_file = save_Gtiff(tvdi, evi_f, zona, opc)
    #######################################

    return tvdi_file

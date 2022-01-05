import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from collections import Counter

parser = argparse.ArgumentParser(description="""    Final de Introducción a Python.
                                                    Agustin Parrotta.
                                                    Analisis criptografico de texto.
                                                    Este script cifra o descifra un mensaje mediante el metodo de Cifrado Afin.
                                                    Para cifrar un mensaje, se deben introducir las claves 'a' y 'b'.
                                                    Para descifrar un mensaje, no es necesario introducir 'a' y 'b'.
                                                    Tambien resuelve el ejercicio extra de descifrar mensaje_cifrado_10
                                                    Tambien obtiene los graficos de distribucion de caracteres para los idiomas Español, Ingles, Aleman y Finlandes""")
parser.add_argument('modo', type=str, help="cifrar, descifrar, extra, idiomas" , choices=['cifrar','descifrar','extra', 'idiomas'])
parser.add_argument("archivo", type=str,help="Nombre del archivo donde se encuentra el mensaje para cifrar o descifrar", nargs='?')
parser.add_argument('a', type=int, nargs='?')
parser.add_argument('b', type=int, nargs='?')
args = parser.parse_args()

def leerArchivo (file):
    """Funcion que lee el archivo file
        Input:
        file: string, direccion del archivo a leer
        Output:
        mensaje: string"""

    with open(file, 'r') as f:
        mensaje = f.read()
    return mensaje

def sonCoprimos (a, b):
    """Funcion que define si a y b son primos entre sí
        Input:
        a, b : integer
        Output:
        bool True si son coprimos, False si no lo son"""

    while b != 0:
        a, b = b, a % b
    if (a == 1 or a == -1):
        return True
    else:
        return False

def moduloInverso(a, n):
    """Funcion que calcula el inverso multiplicativo de a en modulo n
        Input:
        a: integer, numero a calcular el moduloInverso
        n: integer, modulo
        Output:
        integer si existe el valor, None si no existe"""
    # La teoría fue obtenida de: http://www.criptored.upm.es/crypt4you/temas/criptografiaclasica/leccion8.html#apartado2-1
    # Algoritmo: https://www.youtube.com/watch?v=D289EF58Yrw&t=192s

    if not sonCoprimos(a, n):
        return None

    u0, u1, g0 = 1, 0, n
    v0, v1, g1 = 0, 1, a
    while g1 != 0:
        y = g0 // g1
        g1, g0, u1, u0, v1, v0 = (g0 - y * g1), g1, (u0 - y * u1), u1, (v0 - y * v1), v1
    return v0 % n

def coefPearson (x, y):
    """Funcion que calcula el coeficiente de Pearson entre dos vectores
        Input:
        x, y: arrays
        Output:
        coeficiente de Pearson"""

    return stats.pearsonr(x, y)

def diccionario ():
    """Funcion que devuelve una lista donde sus elementos son los 27 caracteres
        Output:
        dicc: ['A', 'B', ..., 'Z', ' ']"""

    dicc = [0] * 27

    for i in range(26):
        dicc[i] = chr(ord('A') + i)

    dicc[26] = ' '

    return dicc

def cambiarEspeciales (s):
    """Funcion que reemplaza los caracteres con signo por sus respectivos sin signos
        Input:
        string, mensaje que puede contener caracteres con signo
        Output:
        s: mensaje sin caracteres con signo"""

    lista1 = np.array(['Á', 'É', 'Í', 'Ó', 'Ú', 'Ä', 'Ë', 'Ï', 'Ö', 'Ü', 'Ñ'])
    lista2 = np.array(['A', 'E', 'I', 'O', 'U', 'A', 'E', 'I', 'O', 'U', 'N'])

    for i in range(len(lista1)):
        s = s.replace(lista1[i], lista2[i])

    return s

def frecuencia (s):
    """Funcion que devuelve un vector con la frecuencia de cada caracter en el mensaje s
        Input:
        s: string, mensaje de entrada
        Output:
        array, donde cada elemento i es el porcentaje de aparición del caracter i"""
    s = s.upper()
    s = cambiarEspeciales(s)

    dicc = diccionario()
    frec = np.empty(27)

    for index, value in enumerate(dicc):
        frec[index] = s.count(value)

    return frec * 100 / np.sum(frec)

def obtenerIdiomas (*fidi):
    """Funcion que devuelve un diccionario con la distribucion de caracteres en distintos idiomas
        Input:
        fidi: tuplas (string idioma, string file). Cada tupla contiene en su primer elemento un string que indica el idioma, y
        en el segundo elemento el string con la direccion de un mensaje en el idioma correspondiente
        con el que se pueda calcular la distribucion
        Output:
        dicc: diccionario. la llave es el string idioma. dicc[idioma] contiene un array con la distribucion de frecuencia"""
    dicc = {}
    for i in range(len(fidi)):
        s = leerArchivo(fidi[i][1])
        vectorFrec = frecuencia(s)
        graficoDistribucion (vectorFrec, fidi[i][1])
        dicc[fidi[i][0]] = vectorFrec
    return dicc

def ordenarLista (lista1, lista2):
    """Funcion que ordena la lista1 de mayor a menor segun la lista2
        Input:
        lista1: lista a ordenar
        lista2: lista segun la cual ordenar
        Output:
        lista: lista ordenada, de mayor a menor, segun la lista2"""

    lista = [x for _,x in sorted(zip(lista2, lista1), reverse = True)]

    return lista

def graficoDistribucion (frec, nombre):
    """Funcion que realiza un grafico de barras con la distribucion de los caracteres
        Input:
        frec: array que indica la frecuencia de aparicion de cada caracter
        nombre: string, nombre del archivo donde se guarda el gráfico
        Output:
        guarda el gráfico realizado en 'nombre.png'
        retorna None"""

    with plt.style.context('default'):
        figure = plt.figure()
        plt.bar(diccionario(), frec, label = nombre.strip('data/').strip('.txt'))
        plt.legend(loc = 'upper left')
        plt.xlabel('Caracteres')
        plt.ylabel('Frecuencia de aparición (%)')
        plt.ylim(0, 20)
        plt.savefig('{}.png'.format(nombre.strip('.txt')), dpi = 200)
        plt.close(figure)

def graficoComparacion (frecMsj, frecIdioma, nombre, idioma):
    """Funcion que realiza un grafico de barras comparando la frecuencia de aparicion
    de caracteres en un idioma y en un mensaje
        Input:
        frecMsj: array que indica la frecuencia de aparicion de cada caracter en un mensaje
        frecIdioma: array que indica la frecuencia de aparicion de cada caracter en un idioma
        idioma: string que indica el idioma en cuestion
        nombre: string, nombre del archivo donde se guarda el gráfico
        Output:
        guarda el gráfico realizado en 'nombre.png'
        retorna None"""

    letrasOrdenadasMsj = ordenarLista(diccionario(), frecMsj)
    letrasOrdenadasIdioma = ordenarLista(diccionario(), frecIdioma)

    frecMsj = -np.sort(-frecMsj)
    frecIdioma = -np.sort(-frecIdioma)

    with plt.style.context('default'):
        figure = plt.figure()
        axis1 = figure.add_subplot(111, axisbelow = True)
        plt.grid()
        plt.ylabel('Frecuencia de aparición (%)')

        axis2 = axis1.twiny()

        axis2.bar(letrasOrdenadasMsj, -frecIdioma, color = 'orange', label = 'Idioma {}'.format(idioma))
        axis1.bar(letrasOrdenadasIdioma, frecMsj,label = nombre)

        axis1.legend(loc = 'upper right')
        axis2.legend(loc = 'lower right')
        axis1.axhline(0, color = 'black')

        axis2.set_yticklabels([str(abs(x)) for x in axis2.get_yticks()])

        plt.xlabel('Caracteres')

        plt.savefig('{}.png'.format(nombre), dpi = 200)
        plt.close(figure)

    return

def graficoCorrelacion (frecCifrado, frecDescifrado, frecIdioma, nombre, idioma):
    """Funcion que realiza un grafico de correlacion entre la frecuencia de aparicion
    de los caracteres en un idioma y en el mensaje
        Input:
        frecCifrado: array que indica la frecuencia de aparicion de cada caracter en el mensaje cifrado
        frecDescifrado: array que indica la frecuencia de aparicion de cada caracter en el mensaje descifrado
        frecIdioma: array que indica la frecuencia de aparicion de cada caracter en un idioma
        idioma: string que indica el idioma en cuestion
        nombre: string, nombre del archivo donde se guarda el gráfico
        Output:
        guarda el gráfico realizado en 'nombre.png'
        retorna None"""

    with plt.style.context('default'):
        figure = plt.figure()

        axis1 = figure.add_subplot(111)
        plt.plot(frecIdioma,frecCifrado, 'o',label = 'Cifrado')
        plt.plot(frecIdioma,frecDescifrado, 'o',color = 'black', label = 'Descifrado')

        c = max(*axis1.get_xlim(),*axis1.get_xlim())
        axis1.set_xlim(-1, c +1)
        axis1.set_ylim(-1, c + 1)
        x = np.linspace(*axis1.get_xlim())
        axis1.plot(x - 1, x - 1, 'k--', label = 'Correlación perfecta')

        plt.xlabel('Frecuencia de aparición en idioma {} (%)'.format(idioma))
        plt.ylabel('Frecuencia de aparición en mensaje (%)')
        plt.legend(loc = 'upper left')
        plt.savefig('{}.png'.format(nombre), dpi = 200)
        plt.close(figure)

    return

def descifradoAnalitico (mensaje, a, b):
    """Funcion que descifra el mensaje analíticamente con las claves a y b
        Input:
        a: integer
        b: integer
        mensaje: string
        Output:
        traduccion: string, mensaje descifrado"""

    mensaje = mensaje.upper()

    traduccion = ''

    diccionario = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '

    for letra in mensaje:
        if diccionario.find(letra) != -1:
            if letra != ' ':
                num = ord(letra) - ord('A')
            else:
                num = 26

            num = ((num - b) * moduloInverso(a, 27)) % 27

            if num != 26:
                traduccion += chr(num + ord('A'))
            else:
                traduccion += ' '

        else:
            traduccion += letra

    return traduccion

def cifrado (mensaje, a, b):
    """Funcion que cifra el mensaje analíticamente con las claves a y b
        Input:
        a: integer
        b: integer
        mensaje: string
        Output:
        traduccion: string, mensaje cifrado"""

    mensaje = mensaje.upper()

    if not sonCoprimos (a, 27):
        print("Error. 'a' debe ser coprimo con 27.")
        return None

    else:
        traduccion = ''

        diccionario = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '

        for letra in mensaje:
            if diccionario.find(letra) != -1:
                if letra != ' ':
                    num = ord(letra) - ord('A')
                else:
                    num = 26


                num = (a * num + b) % 27

                if num != 26:
                    traduccion += chr(num + ord('A'))
                else:
                    traduccion += ' '

            else:
                traduccion += letra

        return traduccion


def desencriptar (file, a = None, b = None, idiomas = None):
    """Funcion que desencripta un mensaje. Si se especifican las claves 'a' y 'b',
    descifra el mensaje con estas claves. Sino, intenta desencriptar comparando con
    la distribucion de caracteres en distintos idiomas
        Input:
        file: string, direccion del archivo que contiene el mensaje a desencriptar
        a: integer
        b: integer
        idiomas: diccionario que contiene las frecuencias de aparicion de caracteres en distintos idiomas
        Output:
        imprime en pantalla las claves 'a' y 'b' obtenidas y el mensaje traducido
        Si no se proporcionaron las claves 'a' y 'b', también imprime el idioma detectado
        y el coeficiente de correlación de Pearson
        retorna None"""
    mensaje = leerArchivo(file)

    if (a != None and b != None):
        if not sonCoprimos (a, 27):
            print("Error. 'a' debe ser coprimo con 27.")
            return None
        else:
            print("Se utilizan las claves a: {} y b: {}".format(a, b))
            print("El mensaje descifrado es:\n", descifradoAnalitico(mensaje, a, b))
            return

    else:

        frecEnc = frecuencia(mensaje)
        dicc = np.arange(27)
        diccEnc = ordenarLista(dicc, frecEnc)

        coefs = [0]*len(idiomas)

        maximo = [0]*6

        for i, j, k, l in ((i1, j1, k1, l1) for i1 in range(4) for j1 in range(4) for k1 in range(4) for l1 in range(4)):
            if (i == k and j == l):
                continue

            for index, value in enumerate(idiomas):
                diccIdioma = ordenarLista(dicc, idiomas[value])

                if not sonCoprimos((diccIdioma[j] - diccIdioma[l]), 27):
                    # print("Error. No se puede calcular el inverso")
                    coefs[index] = (value, 0, 0, 0, 0, 0)
                else:
                    a = ((diccEnc[i] - diccEnc[k]) * moduloInverso((diccIdioma[j] - diccIdioma[l]), 27)) % 27
                    if not sonCoprimos(a, 27):
                        # print("Error. 'a' no tiene inverso")
                        coefs[index] = (value, 0, 0, 0, 0, 0)
                    else:
                        b = (diccEnc[i] - a * diccIdioma[j]) % 27
                        traduccion = descifradoAnalitico(mensaje, a, b)
                        frecTrad = frecuencia(traduccion)
                        aux, _ = coefPearson(frecTrad, idiomas[value])
                        coefs[index] = (value, aux, a, b, traduccion, frecTrad)


            maximo = max(coefs, key=lambda x:x[1])

            if(maximo[1] > 0.8):
                break

        print('El idioma es: ', maximo[0])
        print('El coeficiente de Pearson es: ', maximo[1])
        print('Las claves son: a = {} y b = {}'.format(maximo[2], maximo[3]))
        print('El mensaje descifrado es:\n', maximo[4])

        graficoComparacion(frecEnc, idiomas[maximo[0]], 'Mensaje Cifrado', maximo[0])
        graficoCorrelacion(frecEnc, maximo[5], idiomas[maximo[0]], 'Correlacion', maximo[0])
        return

def frecuenciaPalabras (s):
    """Funcion que encuentra las 20 palabras mas comunes dento del mensaje s
        Input:
        s: string, mensaje
        Output:
        tupla de tuplas (string, integer). La tupla contiene 20 palabras mas comunes.
        Cada elemento es una tupla con la palabra y su frecuencia de aparición"""

    s = s.lower()

    sSplit = s.split()
    x = Counter(sSplit)

    return x.most_common(20)

def probarClaves (maximo, s, masComunesIdioma, x1, x2, y1, y2):

    coef = 0

    for a, b in ((a1, b1) for a1 in range(27) for b1 in range(27)):

        if not sonCoprimos(a, 27):
            continue

        c1 = (a * x1 + b) % 27
        c2 = (a * x2 + b) % 27

        if ((c1 == y1 and c2 == y2) or (c1 == y2 and c2 == y1)):

            traduccion = descifradoAnalitico(s, a, b)

            masComunes = frecuenciaPalabras(traduccion)
            masComunes = [x[0] for x in masComunes]

            contador = 0
            for item in masComunes:
                if item in masComunesIdioma:
                    contador += 1

            if contador > maximo[2]:
                maximo = (a, b, contador, traduccion)

    return maximo


def ejercicioExtra1 (file, idioma, fesp):

    mensaje = leerArchivo(file)
    frecEnc = frecuencia(mensaje)

    dicc = np.arange(27)
    diccEnc = ordenarLista(dicc, frecEnc)
    diccIdioma = ordenarLista(dicc, idioma)

    # Al hacer la distribucion de caracteres del mensaje encriptado, observe que hay dos caracteres muy frecuentes en comparacion al resto. Por esa razon,
    # supuse que era un trabalengua, por lo cual esos dos caracteres tendrian que corresponderse a un espacio y a una vocal

    vocales = [0, 4, 8, 14, 20]
    coef = [0] * len(vocales)
    maximo = [0] * 4

    mensajeIdioma = leerArchivo(fesp)
    masComunesEsp = frecuenciaPalabras(mensajeIdioma)
    masComunesEsp = [x[0] for x in masComunesEsp]

    for i, vocal in enumerate(vocales):
        maximo = probarClaves(maximo, mensaje, masComunesEsp, vocal, diccIdioma[0], diccEnc[0], diccEnc[1])
    print('Las claves son: a={} b={}'.format(maximo[0], maximo[1]))
    print('El mensaje es: ', maximo[3])
    return



fesp = 'quijote_es.txt'
feng = 'quijote_en.txt'
fdeu = 'quijote_de.txt'
ffin = 'quijote_fi.txt'

idiomas = obtenerIdiomas(('Español', fesp), ('Ingles', feng), ('Aleman', fdeu), ('Finlandes', ffin))

if args.modo=="cifrar":
    msjCifrado = cifrado(leerArchivo(args.archivo), args.a, args.b)
    print("El mensaje cifrado con las claves {} y {} es:\n{}".format(args.a, args.b, msjCifrado))
    msjDescifrado = descifradoAnalitico(msjCifrado, args.a, args.b)
    print("El mensaje descifrado con las mismas claves {} y {} es:\n{}".format(args.a, args.b, msjDescifrado))
if args.modo=="descifrar":
    desencriptar(args.archivo, args.a,args.b, idiomas)
if args.modo=="extra":
    ejercicioExtra1('mensaje_cifrado_10.txt', idiomas['Español'], 'quijote_es.txt')
if args.modo=="idiomas":
    obtenerIdiomas(('Español', fesp), ('Ingles', feng), ('Aleman', fdeu), ('Finlandes', ffin))

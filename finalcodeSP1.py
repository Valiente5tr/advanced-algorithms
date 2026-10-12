#Configuracion general
FRAMESHIFT = 13467 # Posicion base 0 donde ocurre el frameshift de ORF1ab
#posicion biologica 13468


#Funciones generales
def leerFASTA(lecturearchive):
    #Complejidad de tiempo: O(n)
    # Se recorren una sola vez todos los caracteres del archivo.
    #Complejidad de espacio: O(n)
    #Se almacena la secuencia completa en memoria.
    # #leerFASTA: funciona para leer el archivo, eliminar el encabezado y regresar la secuencia unida.

    with open(lecturearchive, "r") as archivo: #se cierra el archivo despues de lectura con with
        encabezado = archivo.readline().strip() #para eliminar el encabezado
        lineasSecuencia = [] #lista para lineas de ADN
        #se opta por una lista que concantenar strings varias veces}
        #mucho mas sencillo y eficiente

        #recorrido de todo el archivo
        for linea in archivo:

            linea = linea.strip() #strip para quitar salto de linea o espacio
            #en caso de alguna linea vacia 
            if linea != "":
                lineasSecuencia.append(linea) #se limpia y se añade a lista

        secuencia = "".join(lineasSecuencia) #se unen todas las lineas para dejar una sola cadena
        # de ["ATGC", "TTAC", "GGAA"] a "ATGCTTACGGAA"

    return secuencia
def leerProteinas(proteinArchive):
    #Misma complejidad que leerFASTA: Solamente se recorre la seccion una vez para todos los caracteres del archivo
    proteinas = [] #lista para proteinas econtradas
    #se almacena temporalmente la proteina que se esta leyendo
    nombre = ""
    secuencia = []

    with open(proteinArchive, "r") as archivo:
        for linea in archivo:
            linea = linea.strip()

            if linea == "":
                continue

            if linea.startswith(">"): #> entonces es el encabezado de otra proteina
                if nombre != "": #se guarda la proteina anterior antes de leer la que sigue
                    proteinas.append({"nombre": nombre,"secuencia": "".join(secuencia)})

                nombre = linea[1:] #> eliminamos el caracter y se gurdar el identificador
                secuencia = [] #otra secuencia nueva

            else:
                secuencia.append(linea) #en caso de que no haya > entonces sigue siendo parte de la secuencia

    #se guarda la ultima proteina despues de terminar el archivo (despues no habra un > para finalizarla)
    if nombre != "":
        proteinas.append({"nombre": nombre,"secuencia": "".join(secuencia)})

    return proteinas
def crearTablaCodones():

    #Complejidad de tiempo: O(1) La tabla siempre contiene exactamente los 64 codones posibles.
    # Complejidad de espacio: O(1) El numero de codones posibles es fijo.
    # Consultar un codon en el diccionario: O(1) en promedio.
    # Para esta parte consideramos que varios tipos de codones producen aminoacidos iguales, 
    # entonces cuando aparezca cierta letra no podemos considerar exclusivamente que vino de un codon en especifico 
    # y para eso tenemos que conservar el codon exacto que se encuentra en el genoma


    # Cada llave es un codon formado por 3 bases de ADN.
    # El valor es el aminoacido que produce.
    # el caracter * signifca que es un codon de paro

    tabla = {
        # Fenilalanina (F) y Leucina (L)
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        # Serina (S)
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        # Tirosina (Y) y codones STOP
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        # Cisteina (C), STOP y Triptofano (W)
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",

        # Leucina (L)
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        # Prolina (P)
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        # Histidina (H) y Glutamina (Q)
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        # Arginina (R)
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        # Isoleucina (I) y Metionina (M)
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        # Treonina (T)
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        # Asparagina (N) y Lisina (K)
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        # Serina (S) y Arginina (R)
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        # Valina (V)
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        # Alanina (A)
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        # Acido aspartico (D) y Acido glutamico (E)
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        # Glicina (G)
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G"
    }
    return tabla

#Funciones para PALINDROMOS
def preprocess(cadena):
    resultado = "%"
    
    for caracter in cadena:
        resultado += "#"
        resultado += caracter
    
    resultado += "#$"
    return resultado
def check(char1, char2):
    #PODEMOS UTILIZAR UN DICCIONARIO
    #Comprobacion si dos caracteres son complementarios
    #A-T  y G-C pares en ADN
    match char1:
        case "A":
            if char2 == "T":
                return True
            else:
                return False
        case "T":
            if char2 == "A":
                return True
            else:
                return False
        case "G":
            if char2 == "C":
                return True
            else:
                return False
        case "C":
            if char2 == "G":
                return True
            else:
                return False
        case "#":
            return char2 == "#"
        case _:
            return False
def manacher(cadena):
    #Complejidad de tiempo: O(n)
    #Complejidad de espacio: O(n)
    #El arreglo p almacena el radio del palindromo encontrado alrededor de cada posicion.
    cadena_t = preprocess(cadena)
    n = len(cadena_t)

    p = [0 for i in range(n)]

    center = 0
    right = 0

    maxCenter = 0
    maxLength = 0

    for i in range(1, n - 1):
        if cadena_t[i] != "#":
            continue
            
        mirror = 2 * center - i

        if i < right:
            p[i] = min(right - i, p[mirror])

        while check(cadena_t[i + (1 + p[i])], cadena_t[i - (1 + p[i])]):
            p[i] += 1

        if i + p[i] > right:
            center = i
            right = i + p[i]

        if p[i] > maxLength:
            maxLength = p[i]
            maxCenter = i

    start = (maxCenter - maxLength) // 2

    return cadena[start : start + maxLength]

#Funciones pde PROTEINAS
def traducirMarco(secuencia, desplazamiento, tablaCodones):

    #Complejidad: O(n) Recorremos la secuencia una sola vez y avanzando de 3 en 3.
    aminoacidos = "" #Secuencia final de aminoacidos
    informacionCodones = [] #codon traducido para recuperar posicion en el genoma al encontrar una proteina

    for i in range(desplazamiento, len(secuencia) - 2, 3): #-2 al limite para que siempre existan 3 caracteres

        codon = secuencia[i:i + 3] #se toman 3 nucleotidos
        aminoacido = tablaCodones[codon] #que aminoacido lo produce
        aminoacidos += aminoacido

        # Guardamos info para regresar despues el ADN original
        informacionCodones.append({"codon": codon, "aminoacido": aminoacido,"inicio": i,"fin": i + 2})

    return aminoacidos, informacionCodones
def complementoInverso(secuencia):
    # Complejidad de tiempo: O(n)
    # Complejidad de espacio: O(n)

    complemento = {"A": "T","T": "A","C": "G","G": "C"} #Bases complementarias para el ADN
    secuenciaComplementaria = [] #construccion complementaria

    for base in secuencia:

        nuevaBase = complemento[base] #se busca la base
        secuenciaComplementaria.append(nuevaBase) #se agrega a la nueva

    secuenciaComplementaria = "".join(secuenciaComplementaria) #la lista la volvemos a string
    secuenciaInversa = secuenciaComplementaria[::-1] #se invierte la cadena

    return secuenciaInversa
def construirLPS(patron):
    #Complejidad de tiempo: O(m) Se recorre el patron sin regresar completamente hacia atras.
    #Complejidad de espacio: O(m) Se crea un arreglo LPS del mismo tamaño que el patron.
    #LPS con misma longitud al patron: por posicion se guardan cuantos caracteres se pueden utilizar en caso de diferencia en la busqueda
    lps = [0] * len(patron)
    longitud = 0 #tamanio del prefijo que coincide con un sufijo
    i = 1 #lps[0] #siempre vale 0

    while i < len(patron):

        # Si los caracteres coinciden se aumenta la longitud del prefijo valido
        if patron[i] == patron[longitud]:
            longitud += 1
            lps[i] = longitud
            i += 1

        else:
            if longitud != 0: #en caso de no haber coincidencia pero ya exisitia una previa

                longitud = lps[longitud - 1] #se reutiliza el prefijo mas chico

            else:
                lps[i] = 0    #LPS igual a 0 si tampoco hay coincidencia anterior
                i += 1

    return lps
def KMP(texto, patron):

    if len(patron) == 0: #no se realiza busqueda con patron vacio
        return -1
    lps = construirLPS(patron) #LPS del patron
    i = 0 #para recorrer el texto
    j = 0 #para el patron

    while i < len(texto):
        if texto[i] == patron[j]:
            i += 1  #aumentamos si las cadenas son iguales
            j += 1

            if j == len(patron): #se encontro toda la proteina cuando se llega al final del patron
                return i - j #se regresa la posicion donde inicia coincidencia

        else:
            if j != 0: #ya se habia econtrado parte del patron
                j = lps[j - 1]

            else:
                i += 1 #sin coincidencia
    return -1
def traducirORF1abConFrameshift(genoma, tablaCodones):

    # Las coordenadas de ORF1ab en el archivo biologico
    inicioORF = 265
    posicionFrameshift = FRAMESHIFT
    finORF = 21554

    #construccion manual
    primeraParte = genoma[inicioORF : posicionFrameshift + 1]
    segundaParte = genoma[posicionFrameshift : finORF + 1]
    secuenciaORF = primeraParte + segundaParte
    indicesGenoma = []

    # Indices correspondientes a la primera parte.
    for i in range(inicioORF, posicionFrameshift + 1):
        indicesGenoma.append(i)

    # Indices correspondientes a la segunda parte.
    for i in range(posicionFrameshift, finORF + 1):
        indicesGenoma.append(i)

    aminoacidos = "" #para traduccion final
    informacionCodones = [] #codon generado


    # Traducimos de 3 en 3.
    for i in range(0, len(secuenciaORF) - 2, 3):
        codon = secuenciaORF[i:i + 3]
        aminoacido = tablaCodones[codon]
        aminoacidos += aminoacido

        # Obtenemos las posiciones reales del genoma que forma el codon
        posiciones = [indicesGenoma[i],indicesGenoma[i + 1],indicesGenoma[i + 2]]
        informacionCodones.append({"codon": codon,"aminoacido": aminoacido,"posiciones": posiciones})

        #Justificación: ORF1ab requiere esta función específica debido a un cambio de marco 
        #de lectura ribosomal -1 (-1 ribosomal frameshift) en la posición 13468.
        #Al participar dos veces dicho nucleótido (join(266..13468, 13468..21555)), 
        #se altera la agrupación de los codones. Para garantizar máxima precisión y evitar errores de traducción, 
        #construimos el gen forzando las coordenadas anotadas exactas del genoma en lugar de buscar 
        #la secuencia deslizante (TTTAAAC) de forma automatizada.

        #Complejidad de tiempo: O(n) Se construye la secuencia
        #y se recorre Complejidad de espacio: O(n) se guarda la secuencia reconstruida, 
        #los indices originales, la traduccion, la informacion de los codones

    return aminoacidos, informacionCodones
def generarTraducciones(genoma, tablaCodones):
    traducciones = [] #se guardan las traducciones para buscar despues las proteinas

    #hebra original
    for frame in range(3):
        aminoacidos, informacion = traducirMarco(genoma,frame,tablaCodones)
        #se normalliza la informacion, la traduccion tendra lista para posiciones
        codones = []

        for dato in informacion:
            codones.append({"codon": dato["codon"],"aminoacido": dato["aminoacido"],"posiciones": [dato["inicio"],dato["inicio"] + 1, dato["fin"]]})
        traducciones.append({"nombre": "Frame +" + str(frame),"hebra": "directa", "frame": frame,"aminoacidos": aminoacidos,"codones": codones}) #hebra directa es original

    #Complemento inverso se requiere:
    inversa = complementoInverso(genoma)
    longitudGenoma = len(genoma)

    #hebra inversa
    for frame in range(3):
        aminoacidos, informacion = traducirMarco(inversa,frame,tablaCodones)
        codones = []

        for dato in informacion:

            #los indices de la inversa no son iguales al del genoma original
            #entonces si el genoma tiene long n entonces posicon original = n-1-posicion inversa

            posicion1 = longitudGenoma - 1 - dato["inicio"]
            posicion2 = longitudGenoma - 1 - (dato["inicio"] + 1)
            posicion3 = longitudGenoma - 1 - dato["fin"]

            codones.append({"codon": dato["codon"],"aminoacido": dato["aminoacido"], "posiciones": [posicion1,posicion2,posicion3]})

        traducciones.append({"nombre": "Frame -" + str(frame),"hebra": "inversa","frame": frame,"aminoacidos": aminoacidos,"codones": codones})


    #CASO PARA EL ORFab
    aminoacidosORF, codonesORF = traducirORF1abConFrameshift(genoma,tablaCodones)

    traducciones.append({"nombre": "ORF1ab con frameshift","hebra": "directa", "frame": "especial", "aminoacidos": aminoacidosORF, "codones": codonesORF})

    #Justificación: Centraliza y unifica la generación de los 7 marcos de lectura: 
    #3 normales, 3 de la hebra complementaria inversa y el caso especial de ORF1ab.
    #Separar esto de buscarProteina() respeta el principio de responsabilidad única: 
    #esta función limpia y mapea las coordenadas (corrigiendo el desfase de índices en la hebra inversa),
    #dejando una lista uniforme lista para que el algoritmo KMP se dedique exclusivamente a la búsqueda.
    return traducciones
def buscarProteina(proteina, traducciones):

    #Complejidad de tiempo: O(n + m) n = cantidad de aminoacidos de las traducciones. m = longitud de la proteina buscada.
    #KMP realiza la busqueda en O(n + m). Aunque se esten revisando varias traducciones,el numero de reading frames es constante.

    #Datos de la proteina a buscar
    nombre = proteina["nombre"]
    secuenciaProteina = proteina["secuencia"]

    #se revisan las traducciones
    for traduccion in traducciones:

        textoAminoacidos = traduccion["aminoacidos"]
        #utilizando el KMP se busca la secuencia completa de la proteina
        posicion = KMP(textoAminoacidos,secuenciaProteina)
        #no esta en la traduccion
        if posicion == -1:
            continue

        # Si llegamos aqui, encontramos la proteina.
        #posicoion es el indice del aminoacido y donde empieza en la traduccion
        indicePrimerAminoacido = posicion
        indiceUltimoAminoacido = (posicion + len(secuenciaProteina) - 1)

        codonesProteina = traduccion["codones"][indicePrimerAminoacido :indiceUltimoAminoacido + 1] #se recuperan los codones
        #posiciones a recuperar

        posicionesGenoma = []

        for datoCodon in codonesProteina:

            for indice in datoCodon["posiciones"]:

                posicionesGenoma.append(indice)

        #limite de la region del genoma , min y max  porque en la inversa las posiciones pueden ser de mayor a menor
        inicioGenoma = min(posicionesGenoma)
        finGenoma = max(posicionesGenoma)


        primerosAminoacidos = secuenciaProteina[:4] #primeros 4 aminoacidos
        primerosCodones = [] #los primeros cuatro codones

        limite = min(4, len(codonesProteina))

        for i in range(limite):
            primerosCodones.append(codonesProteina[i]["codon"])

        #se guardan los resultados
        resultado = {"nombre": nombre,"traduccion": traduccion["nombre"],"hebra": traduccion["hebra"],"inicio": inicioGenoma + 1,"fin": finGenoma + 1,"primerosAminoacidos": primerosAminoacidos,"primerosCodones": primerosCodones}
        #sumamos 1 para las posiciones biologicas

        return resultado


    # Si revisamos todas las traducciones
    # y nunca encontramos la proteina.
    return None
def mostrarResultado(resultado):

    #si es none no se encuentra en las traducciones
    if resultado is None:
        print("Proteina no encontrada")
        print()
        return

    print("Proteina:", resultado["nombre"]) #nombre
    print("Indices en el genoma:",resultado["inicio"],"-",resultado["fin"])
    print("Primeros 4 aminoacidos:",resultado["primerosAminoacidos"])

    print("Codones asociados:")

    aminoacidos = resultado["primerosAminoacidos"]
    codones = resultado["primerosCodones"]

    for i in range(len(codones)):
        print(aminoacidos[i],"->",codones[i])


    #COMPROBACION DE COMO SE ENCONTRO
    print("Reading frame:",resultado["traduccion"])
    print("Hebra:",resultado["hebra"])
    print()

#Funciones para COMPARACION
def computeLPS(pattern, M, lps):
    # Complejidad de tiempo: O(m)
    # Complejidad de espacio adicional: O(1)
    longitud = 0
    lps[0] = 0
    i = 1

    while i < M:
        if pattern[i]==pattern[longitud]:
            longitud += 1
            lps[i] = longitud
            i += 1
        else:
            if longitud != 0:
                longitud = lps[longitud - 1]
            else:
                lps[i] = 0
                i += 1
def searchIndices(pattern, txt):
    # Complejidad de tiempo: O(n + m). Espacio: O(m) del arreglo LPS
    M = len(pattern)
    N = len(txt)

    if M == 0:
        return []

    lps = [0]*M
    j = 0
    indices = []

    computeLPS(pattern, M, lps)

    i = 0
    while i < N:
        if pattern[j] == txt[i]:
            i += 1
            j += 1

        if j == M:
            indices.append((i - j, i - 1))
            j = lps[j - 1]

        elif i < N and pattern[j] != txt[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return indices
def coinciden(a, ia, b, ib, k):
    # Compara k caracteres de a desde ia contra b desde ib, caracter por caracter
    # Complejidad: O(k)
    if ia + k > len(a) or ib + k > len(b):
        return False
    for t in range(k):
        if a[ia + t] != b[ib + t]:
            return False
    return True
def anclar(wuhan, texas, k):
    # Encuentra donde se alinean los inicios de ambos genomas buscando con KMP
    # un k-mer de Texas en Wuhan (salta k-mers con bases ambiguas como N).
    # Si el k-mer de Texas en s esta en Wuhan en p, texas[0] corresponde a wuhan[p - s].
    # Complejidad: O(n + m) por k-mer probado
    for s in range(0, len(texas) - k, k):
        kmer = texas[s:s + k]
        #ignorar regiones con bases ambiguas
        if "N" in kmer:
            continue

        encontrados = searchIndices(texas[s:s + k], wuhan)
        if len(encontrados) > 0:
            desfase = encontrados[0][0] - s
            if desfase >= 0:
                return desfase, 0
            return 0, -desfase
    return None
def reanclar(wuhan, i, texas, j, k, ventana):
    # Tras una diferencia que no es sustitucion, busca con KMP el menor salto que
    # vuelva a alinear: delecion en Texas (su k-mer aparece mas adelante en Wuhan)
    # o insercion en Texas (el k-mer de Wuhan aparece mas adelante en Texas).
    # Regresa (tipo, salto) o None. Complejidad: O(ventana + k)
    mejor = None
    enW = searchIndices(texas[j:j + k], wuhan[i + 1:i + 1 + ventana + k])
    if len(enW) > 0:
        mejor = ("delecion", enW[0][0] + 1)
    enT = searchIndices(wuhan[i:i + k], texas[j + 1:j + 1 + ventana + k])
    if len(enT) > 0:
        salto = enT[0][0] + 1
        if mejor is None or salto < mejor[1]:
            mejor = ("insercion", salto)
    return mejor
def compararGenomas(wuhan, texas, k, ventana):
    # Recorre ambos genomas alineados y regresa la lista de diferencias:
    # {posW, posT, baseW, baseT, extremo}
    # Complejidad: O(n + m) en el recorrido + O(ventana) por cada indel
    # Cuando aparece una posible insercion o delecion se realiza
    # una busqueda adicional dentro de una ventana limitada.
    ancla = anclar(wuhan, texas, k)
    if ancla is None:
        return []
    i, j = ancla
    diferencias = []

    # Lo que sobra al inicio antes del ancla
    if i > 0:
        diferencias.append({"posW": 0, "posT": 0, "baseW": wuhan[0:i], "baseT": "-", "extremo": True})
    if j > 0:
        diferencias.append({"posW": 0, "posT": 0, "baseW": "-", "baseT": texas[0:j], "extremo": True})

    while i < len(wuhan) and j < len(texas):
        if wuhan[i] == texas[j]:
            i += 1
            j += 1
            continue

        if coinciden(wuhan, i + 1, texas, j + 1, k):
            # Sustitucion: lo que sigue vuelve a coincidir
            diferencias.append({"posW": i, "posT": j, "baseW": wuhan[i], "baseT": texas[j], "extremo": False})
            i += 1
            j += 1
        else:
            r = reanclar(wuhan, i, texas, j, k, ventana)
            if r is None:
                # Varias sustituciones seguidas: se registra esta y se sigue
                diferencias.append({"posW": i, "posT": j, "baseW": wuhan[i], "baseT": texas[j], "extremo": False})
                i += 1
                j += 1
            elif r[0] == "delecion":
                diferencias.append({"posW": i, "posT": j, "baseW": wuhan[i:i + r[1]], "baseT": "-", "extremo": False})
                i += r[1]
            else:
                diferencias.append({"posW": i, "posT": j, "baseW": "-", "baseT": texas[j:j + r[1]], "extremo": False})
                j += r[1]

    # Lo que sobra al final
    if i < len(wuhan):
        diferencias.append({"posW": i, "posT": j, "baseW": wuhan[i:], "baseT": "-", "extremo": True})
    if j < len(texas):
        diferencias.append({"posW": i, "posT": j, "baseW": "-", "baseT": texas[j:], "extremo": True})

    return diferencias
def rangosProteinas(proteinas, traducciones):
    # Usa el resultado de buscarProteina para obtener el rango de cada proteina en el genoma (hebra directa, marcos normales; 
    # ORF1AB con frameshift se maneja con genes). Regresa (inicio, fin) en base 0.
    # Complejidad: O(n + m) cada una
    rangos = []
    for proteina in proteinas:
        r = buscarProteina(proteina, traducciones)
        if r is not None and r["hebra"] == "directa" and r["traduccion"].startswith("Frame"):
            rangos.append((r["inicio"] - 1, r["fin"] - 1))  # buscarProteina reporta en base 1
    return rangos
def inicioRegion(pos, genes, rangos):
    # Inicio (base 0) de la region codificante que contiene pos: primero los genes M, S y ORF1AB, si no, las proteinas de la parte 3. None si no codifica.
    # Complejidad: O(g + p)
    for nombre in genes:
        for (inicio, fin) in genes[nombre]:
            if inicio <= pos <= fin:
                if nombre == "ORF1AB" and pos > FRAMESHIFT:
                    # despues del frameshift el marco reinicia aqui
                    return FRAMESHIFT
                return inicio
    for (inicio, fin) in rangos:
        if inicio <= pos <= fin:
            return inicio
    return None
def analizarDiferencia(d, wuhan, texas, genes, rangos, tabla):
    # Agrega el codon afectado y el aminoacido de cada version.
    # El codon empieza en inicio + 3 * ((pos - inicio) // 3), respetando el marco de lectura.
    # Complejidad: O(g + p)
    d["codonW"] = None
    if d["extremo"]:
        return d
    inicio = inicioRegion(d["posW"], genes, rangos)
    if inicio is None:
        return d
    ic = inicio + 3 * ((d["posW"] - inicio) // 3)
    inicioT = d["posT"] - (d["posW"] - ic)   # mismo codon en Texas

    codonW = wuhan[ic:ic + 3]
    codonT = texas[inicioT:inicioT + 3]
    if len(codonW) != 3 or len(codonT) != 3: #comprobacion de que se encuentran completos
            return d

    d["codonW"] = wuhan[ic:ic + 3]
    d["codonT"] = texas[inicioT:inicioT + 3]
    d["aaW"] = tabla.get(d["codonW"], "X")   # X: codon con bases ambiguas (N, etc.)
    d["aaT"] = tabla.get(d["codonT"], "X")
    return d
def mostrarDiferencias(diferencias):
    # Imprime indice, codones afectados y aminoacido de cada version. Complejidad: O(d)
    print("Diferencias entre Wuhan 2019 y Texas 2020")
    print("Total de diferencias:", len(diferencias))
    print()
    for d in diferencias:
        print("Indice", d["posW"] + 1, "(Wuhan) /", d["posT"] + 1, "(Texas):", d["baseW"], "->", d["baseT"])
        if d["codonW"] is None:
            print("  Fuera de una region codificante: no produce aminoacido")
        else:
            if d["aaW"] == d["aaT"]:
                resultado = "mismo aminoacido"
            else:
                resultado = "aminoacido diferente"
            print("  Codon Wuhan:", d["codonW"], "(" + d["aaW"] + ") | Codon Texas:", d["codonT"], "(" + d["aaT"] + ") ->", resultado)
        print()
def cadenaMasLargaAlineada(wuhan,texas,diferencias):
    #Se busca la cadena consecutiva mas larga sin diferencias considerando que el alineamiento ya fue calculado
    #COMPLEJIDAD o(n) considerando n como la cantidad de diferencias encontradas
    cursorW = 0
    cursorT = 0
    mejorLongitud = 0
    mejorInicioW = 0
    mejorInicioT = 0

    for diferencia in diferencias:
        posW = diferencia["posW"]
        posT = diferencia["posT"]
        longitudW = posW - cursorW #misma candtidad de bases iguales desde la ultima diferencia
        longitudT = posT - cursorT
        longitudActual = min(longitudW, longitudT) #como estan alineados el tramo debe de estar en ambos

        if longitudActual > mejorLongitud: #en caso de cadena mayor guardamos su posicion
            mejorLongitud = longitudActual
            mejorInicioW = cursorW
            mejorInicioT = cursorT

        #Insercion en Texas:
        #Wuhan no consume ninguna base.
        if diferencia["baseW"] == "-":
            cursorW = posW
            cursorT = (posT + len(diferencia["baseT"]))

        # Delecion en Texas y avanzamos las bases correspondientes en Wuhan.
        elif diferencia["baseT"] == "-":
            cursorW = (posW + len(diferencia["baseW"]))
            cursorT = posT
        else: #ambos genomas toman bases
            cursorW = (posW+ len(diferencia["baseW"]))
            cursorT = (posT+ len(diferencia["baseT"]))

    #Despues de la ultima diferencia todavia puedeexistir la cadena sin cambios mas larga.
    longitudW = len(wuhan) - cursorW
    longitudT = len(texas) - cursorT
    longitudActual = min(longitudW, longitudT)

    if longitudActual > mejorLongitud:
        mejorLongitud = longitudActual
        mejorInicioW = cursorW
        mejorInicioT = cursorT

    return (mejorLongitud,mejorInicioW,mejorInicioT)
#Primer caso para Encontrar cadena mas larga (Desechado)
def subcadenaComunMasLarga(s1, s2):
    #SOLUCION GENERAL
    # Longest common substring con programacion dinamica.
    # Recurrencia: T[i][j] = 1 + T[i-1][j-1] si s1[i-1] == s2[j-1], 0 en otro caso.
    # Complejidad de tiempo: O(n * m)  (1-2 minutos con genomas de 30,000 bases)
    # Complejidad de espacio: O(m) en lugar de O(n * m)
    m = len(s2)
    # fila i-1 (primera fila en 0 como padding)
    prev = [0] * (m + 1)
    mejor = 0
    finS1 = 0
    finS2 = 0
    for i in range(1, len(s1) + 1):
        # primera columna en 0 como padding
        current = [0] * (m + 1)
        caracter = s1[i - 1]
        for j in range(1, m + 1):
            if caracter == s2[j - 1]:
                # extiende la diagonal
                current[j] = prev[j - 1] + 1
                if current[j] > mejor:
                    mejor = current[j]
                    finS1 = i
                    finS2 = j
            # si no coinciden, current[j] se queda en 0
        prev = current
    # inicio = fin - longitud (base 0)
    return mejor, finS1 - mejor, finS2 - mejor
def mostrarCadenaMasLarga(wuhan, texas, diferencias):
    # Imprime la cadena mas larga sin diferencias. Complejidad: O(n * m) por la DP
    L, iw, it = cadenaMasLargaAlineada(wuhan, texas, diferencias)
    print("Cadena mas larga sin diferencias:", L, "bases")
    print("  Wuhan:", iw + 1, "-", iw + L, "| Texas:", it + 1, "-", it + L)
    print("Primero 50 caracteres", wuhan[iw:iw + 50])

#funciones de control de ejecucion
def ejecutarGenes(genes, genoma):

    print("PUNTO 1 - BUSQUEDA DE GENES")
    print("\n")

    for nombre, secuencia in genes.items():
        apariciones = searchIndices(secuencia, genoma)

        for inicio, fin in apariciones:

            print("Gen:", nombre)
            print("Indices:",inicio + 1,"-",fin + 1) #empieza en 0 se suma 1 para las coordenadas
            print("Primeros 12 caracteres: ",secuencia[:12])
            print()
def ejecutarPalindromos(genes):

    print("PUNTO 2 - PALINDROMOS")
    print("\n")

    for nombre, secuencia in genes.items():
        palindromo = manacher(secuencia)

        print("Gen:", nombre)
        print("Longitud:", len(palindromo))
        print("Palindromo:", palindromo)
        print()

        nombreArchivo = "palindromo-" + nombre + ".txt" #palindromo archivo
        #guardarlo para mostrar la longitud

        with open(nombreArchivo, "w") as archivo:
            archivo.write(palindromo)
        print("Se guardo en: ", nombreArchivo)
        print()
def ejecutarProteinas(proteinas, traducciones):

    print("PUNTO 3 - PROTEINAS")
    print("\n")
    encontradas = 0

    for proteina in proteinas:

        resultado = buscarProteina(proteina,traducciones)
        mostrarResultado(resultado)
        if resultado is not None:
            encontradas += 1


    print("Proteinas encontradas:",encontradas,"/",len(proteinas))
def ejecutarComparacion(wuhan,texas,genes,proteinas,traducciones,tablaCodones,k,ventana):
    print("PUNTO 4 - COMPARACION DE GENOMAS")
    print("\n")

    posicionesGenes = {} #genes en Wuhan
    for nombre, secuencia in genes.items():
        posicionesGenes[nombre] = searchIndices(secuencia,wuhan)

    rangos = rangosProteinas(proteinas,traducciones)

    print("Longitud Wuhan:",len(wuhan))
    print("Longitud Texas:",len(texas))
    print()
    diferencias = compararGenomas(wuhan,texas,k,ventana)

    for diferencia in diferencias:
        analizarDiferencia(diferencia,wuhan,texas,posicionesGenes,rangos,tablaCodones)

    mostrarDiferencias(diferencias)
    cadenaMasLargaAlineada(wuhan,texas, diferencias)


def main():

    print("ANALISIS DE GEMONA SARS-CoV-2")

    #lectura y almacenamiento de archivos:
    wuhan = leerFASTA("SARS-COV-2-MN908947.3.txt")
    texas = leerFASTA("SARS-COV-2-MT106054.1.txt")

    genM = leerFASTA("gen-M.txt")
    genS = leerFASTA("gen-S.txt")
    genORF1AB = leerFASTA("gen-ORF1AB.txt")
    proteinas = leerProteinas("seq-proteins.txt")

    #variablesconfiguration
    K = 20 #tam k-mer utilizado para comparacion en los genomas
    VENTANA = 100 #intento para alinear genomas

    #general structures
    genes = {"M": genM, "S":genS, "ORF1AB": genORF1AB}
    tablaCodones = crearTablaCodones()

    #generacion de traducciones para reutilizarlas
    traducciones = generarTraducciones(wuhan, tablaCodones)

    #PUNTO1
    ejecutarGenes(genes,wuhan)
    #PPUNTO2
    ejecutarPalindromos(genes)
    #PUNTO3
    ejecutarProteinas(proteinas,traducciones)
    #PUNTO4
    ejecutarComparacion(wuhan,texas,genes,proteinas,traducciones,tablaCodones,K,VENTANA)


    print("Fin analisis")



#solamente cuando sea el programa principal
if __name__ == "__main__":
    main()

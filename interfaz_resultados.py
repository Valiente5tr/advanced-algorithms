from pathlib import Path
import threading
import tkinter as tk
from tkinter import ttk, messagebox

import finalcodeSP1 as sp

BASE_DIR = Path(__file__).resolve().parent

ARCHIVOS = {
    "wuhan": "SARS-COV-2-MN908947.3.txt",
    "texas": "SARS-COV-2-MT106054.1.txt",
    "M": "gen-M.txt",
    "S": "gen-S.txt",
    "ORF1AB": "gen-ORF1AB.txt",
    "proteinas": "seq-proteins.txt",
}

K = 20
VENTANA = 100

class InterfazResultados:
    def __init__(self, root):
        self.root = root
        self.root.title("Analisis de genoma SARS-CoV-2")
        self.root.geometry("1280x760")
        self.root.minsize(1000, 650)

        self.wuhan = None
        self.texas = None
        self.diferencias = []

        self.crear_interfaz()
        self.root.after(100, self.cargar_resultados)

    def crear_interfaz(self):
        titulo = ttk.Label(self.root,text="Analisis de genoma SARS-CoV-2",font=("Segoe UI", 18, "bold"))
        titulo.pack(pady=(15, 5))
        self.estado = ttk.Label(self.root,text="Cargando resultados...",font=("Segoe UI", 10))
        self.estado.pack(pady=(0, 10))

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=15, pady=(0, 15))

        self.tabGenes = ttk.Frame(self.notebook)
        self.tabPalindromos = ttk.Frame(self.notebook)
        self.tabProteinas = ttk.Frame(self.notebook)
        self.tabComparacion = ttk.Frame(self.notebook)

        self.notebook.add(self.tabGenes, text="1. Genes")
        self.notebook.add(self.tabPalindromos, text="2. Palindromos")
        self.notebook.add(self.tabProteinas, text="3. Proteinas")
        self.notebook.add(self.tabComparacion, text="4. Comparacion")

        self.crear_tabla_genes()
        self.crear_tabla_palindromos()
        self.crear_tabla_proteinas()
        self.crear_tabla_comparacion()

    def crear_treeview(self, padre, columnas, anchos):
        contenedor = ttk.Frame(padre)
        contenedor.pack(fill="both", expand=True, padx=15, pady=15)

        tabla = ttk.Treeview(contenedor, columns=columnas, show="headings")

        for columna in columnas:
            tabla.heading(columna, text=columna)
            tabla.column(columna, width=anchos.get(columna, 120), anchor="center")

        scrollVertical = ttk.Scrollbar(contenedor, orient="vertical", command=tabla.yview)
        scrollHorizontal = ttk.Scrollbar(contenedor, orient="horizontal", command=tabla.xview)

        tabla.configure(yscrollcommand=scrollVertical.set,xscrollcommand=scrollHorizontal.set)

        tabla.grid(row=0, column=0, sticky="nsew")
        scrollVertical.grid(row=0, column=1, sticky="ns")
        scrollHorizontal.grid(row=1, column=0, sticky="ew")

        contenedor.rowconfigure(0, weight=1)
        contenedor.columnconfigure(0, weight=1)

        return tabla

    def crear_tabla_genes(self):
        columnas = ("Gen", "Inicio", "Fin", "Primeros 12 caracteres")
        anchos = {"Gen": 120,"Inicio": 130,"Fin": 130,"Primeros 12 caracteres": 260}
        self.tablaGenes = self.crear_treeview(self.tabGenes, columnas, anchos)

    def crear_tabla_palindromos(self):
        columnas = ("Gen", "Longitud", "Palindromo", "Archivo generado")
        anchos = {"Gen": 110,"Longitud": 110,"Palindromo": 380,"Archivo generado": 250}
        self.tablaPalindromos = self.crear_treeview(self.tabPalindromos, columnas, anchos)

    def crear_tabla_proteinas(self):
        columnas = ("Proteina","Inicio","Fin", "Primeros 4 AA","Codones asociados","Reading frame","Hebra")
        anchos = {"Proteina": 130, "Inicio": 95,"Fin": 95,"Primeros 4 AA": 120,"Codones asociados": 310,"Reading frame": 180,"Hebra": 110}
        self.tablaProteinas = self.crear_treeview(self.tabProteinas, columnas, anchos)

    def crear_tabla_comparacion(self):
        marcoSuperior = ttk.Frame(self.tabComparacion)
        marcoSuperior.pack(fill="x", padx=15, pady=(15, 0))

        self.infoComparacion = ttk.Label(marcoSuperior, text="Esperando datos...")
        self.infoComparacion.pack(side="left")

        self.botonCadena = ttk.Button(marcoSuperior,text="Calcular cadena mas larga",command=self.calcular_cadena_mas_larga, state="disabled")
        self.botonCadena.pack(side="right")

        columnas = ( "Indice Wuhan", "Indice Texas", "Cambio","Frame", "Codon Wuhan", "AA Wuhan", "Codon Texas","AA Texas", "Resultado")
        anchos = { "Indice Wuhan": 110,"Indice Texas": 110,"Cambio": 220,"Frame": 70,"Codon Wuhan": 110, "AA Wuhan": 90,"Codon Texas": 110,  "AA Texas": 90,"Resultado": 180}
        self.tablaComparacion = self.crear_treeview(self.tabComparacion, columnas, anchos)

        self.infoCadena = ttk.Label( self.tabComparacion,text="Cadena mas larga sin diferencias: pendiente",justify="left",font=("Consolas", 10) )
        self.infoCadena.pack(fill="x", padx=15, pady=(0, 15))

    def ruta(self, nombre):
        return str(BASE_DIR / nombre)

    def cargar_resultados(self):
        try:
            self.wuhan = sp.leerFASTA(self.ruta(ARCHIVOS["wuhan"]))
            self.texas = sp.leerFASTA(self.ruta(ARCHIVOS["texas"]))
            genM = sp.leerFASTA(self.ruta(ARCHIVOS["M"]))
            genS = sp.leerFASTA(self.ruta(ARCHIVOS["S"]))
            genORF1AB = sp.leerFASTA(self.ruta(ARCHIVOS["ORF1AB"]))
            proteinas = sp.leerProteinas(self.ruta(ARCHIVOS["proteinas"]))

            genes = {"M": genM, "S": genS, "ORF1AB": genORF1AB}
            tablaCodones = sp.crearTablaCodones()
            traducciones = sp.generarTraducciones(self.wuhan, tablaCodones)

            # PUNTO 1
            for nombre, secuencia in genes.items():
                apariciones = sp.searchIndices(secuencia, self.wuhan)
                for inicio, fin in apariciones:
                    self.tablaGenes.insert( "", "end", values=(nombre, inicio + 1, fin + 1, secuencia[:12]))

            # PUNTO 2
            for nombre, secuencia in genes.items():
                palindromo = sp.manacher(secuencia)
                nombreArchivo = "palindromo-" + nombre + ".txt"
                rutaArchivo = BASE_DIR / nombreArchivo

                with open(rutaArchivo, "w", encoding="utf-8") as archivo:
                    archivo.write(palindromo)

                self.tablaPalindromos.insert( "", "end",values=(nombre, len(palindromo), palindromo, nombreArchivo) )

            # PUNTO 3
            for proteina in proteinas:
                resultado = sp.buscarProteina(proteina, traducciones)

                if resultado is None:
                    self.tablaProteinas.insert("", "end",  values=(proteina["nombre"], "-", "-", "No encontrada", "-", "-", "-"))
                    continue

                codonesTexto = []
                for aminoacido, codon in zip(
                    resultado["primerosAminoacidos"],
                    resultado["primerosCodones"]
                ):
                    codonesTexto.append(aminoacido + " -> " + codon)

                self.tablaProteinas.insert("", "end",
                    values=(
                        resultado["nombre"],
                        resultado["inicio"],
                        resultado["fin"],
                        resultado["primerosAminoacidos"],
                        " | ".join(codonesTexto),
                        resultado["traduccion"],
                        resultado["hebra"]
                    )
                )

            # PUNTO 4
            posicionesGenes = {}
            for nombre, secuencia in genes.items():
                posicionesGenes[nombre] = sp.searchIndices(secuencia, self.wuhan)

            rangos = sp.rangosProteinas(proteinas, traducciones)

            self.diferencias = sp.compararGenomas( self.wuhan,  self.texas,  K, VENTANA)

            for diferencia in self.diferencias:
                esSustitucion = (
                    not diferencia["extremo"]
                    and len(diferencia["baseW"]) == 1
                    and len(diferencia["baseT"]) == 1
                    and diferencia["baseW"] != "-"
                    and diferencia["baseT"] != "-"
                )

                cambio = diferencia["baseW"] + " -> " + diferencia["baseT"]

                if not esSustitucion:

                    self.tablaComparacion.insert(
                        "",
                        "end",
                        values=(
                            diferencia["posW"],
                            diferencia["posT"],
                            cambio,
                            "-",
                            "-",
                            "-",
                            "-",
                            "-",
                            "Indel / extremo"
                        )
                    )

                    continue

                desfase = diferencia["posT"] - diferencia["posW"]

                for frame in range(3):

                    inicioCodonW = diferencia["posW"] - (
                        (diferencia["posW"] - frame) % 3
                    )

                    inicioCodonT = inicioCodonW + desfase

                    if inicioCodonW < 0 or inicioCodonT < 0:
                        continue

                    codonW = self.wuhan[inicioCodonW:inicioCodonW + 3]
                    codonT = self.texas[inicioCodonT:inicioCodonT + 3]

                    if len(codonW) != 3 or len(codonT) != 3:
                        continue

                    aaW = tablaCodones.get(codonW, "X")
                    aaT = tablaCodones.get(codonT, "X")

                    if aaW == aaT:
                        resultadoTexto = "Mismo aminoacido"
                    else:
                        resultadoTexto = "Aminoacido diferente"

                    self.tablaComparacion.insert(
                        "",
                        "end",
                        values=(
                            diferencia["posW"],
                            diferencia["posT"],
                            cambio,
                            frame,
                            codonW,
                            aaW,
                            codonT,
                            aaT,
                            resultadoTexto
                        )
                    )


            self.infoComparacion.config(
                text=(
                    "Wuhan: " + str(len(self.wuhan)) + " bases    |    "
                    "Texas: " + str(len(self.texas)) + " bases    |    "
                    "Diferencias: " + str(len(self.diferencias))
                )
            )

            self.botonCadena.config(state="normal")
            self.estado.config(text="Resultados cargados correctamente.")

        except FileNotFoundError as error:
            self.estado.config(text="No se pudieron cargar los archivos.")
            messagebox.showerror(
                "Archivo no encontrado",
                "Verifica que todos los archivos .txt esten en la misma carpeta.\n\n" + str(error)
            )

        except Exception as error:
            self.estado.config(text="Ocurrio un error.")
            messagebox.showerror("Error", str(error))

    def calcular_cadena_mas_larga(self):
        self.botonCadena.config(state="disabled")
        self.infoCadena.config(text="Calculando cadena mas larga... Esto puede tardar.")

        hilo = threading.Thread(target=self.proceso_cadena_mas_larga,daemon=True)
        hilo.start()

    def proceso_cadena_mas_larga(self):
        try:
            L, inicioW, inicioT = sp.cadenaMasLargaAlineada(self.wuhan, self.texas, self.diferencias)
            primeros50 = self.wuhan[inicioW:inicioW + min(50, L)]

            self.root.after(
                0,
                self.mostrar_resultado_cadena,
                L,
                inicioW,
                inicioT,
                primeros50
            )

        except Exception as error:
            self.root.after(0, self.mostrar_error_cadena, str(error))

    def mostrar_resultado_cadena(self, longitud, inicioW, inicioT, primeros50):
        texto = (
            "Cadena mas larga sin diferencias: " + str(longitud) + " bases\n"
            "Wuhan: " + str(inicioW + 1) + " - " + str(inicioW + longitud)
            + "    |    Texas: " + str(inicioT + 1) + " - " + str(inicioT + longitud)
            + "\nPrimeros 50 caracteres: " + primeros50
        )

        self.infoCadena.config(text=texto)
        self.botonCadena.config(state="normal")

    def mostrar_error_cadena(self, error):
        self.infoCadena.config(text="No fue posible calcular la cadena.")
        self.botonCadena.config(state="normal")
        messagebox.showerror("Error", error)


def main():
    root = tk.Tk()
    InterfazResultados(root)
    root.mainloop()


if __name__ == "__main__":
    main()

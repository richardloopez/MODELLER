#!/bin/bash

# Bucle infinito para contar los archivos y mostrarlos continuamente
while true; do
    # Cuenta los archivos con extensión .pdb en el directorio actual
    # La opción -f en find asegura que solo se cuenten archivos (no directorios)
    # wc -l cuenta las líneas (que corresponden a los archivos encontrados)
    numero_archivos=$(find . -maxdepth 1 -type f -name "*.pdb" | wc -l)

    # Imprime el resultado
    echo "Número de archivos .pdb en el directorio actual: $numero_archivos"

    # Espera 2 segundos antes de repetir el conteo
    sleep 2
done

# NOTA: Este script se ejecutará infinitamente hasta que lo detengas manualmente
# (generalmente presionando Ctrl+C).

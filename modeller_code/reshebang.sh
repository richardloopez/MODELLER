#!/bin/bash
for f in *.py; do
  # Solo si existe el archivo y empieza con shebang
  if [ -f "$f" ]; then
    # Comprueba si la primera línea es shebang (empieza por #!)
    firstline=$(head -n 1 "$f")
    if [[ $firstline == \#\!* ]]; then
      # Reemplaza la primera línea por el shebang correcto
      sed -i '1s|.*|#!/usr/bin/env python3|' "$f"
      echo "Shebang actualizado en $f"
    else
      # Si no hay shebang, lo añade arriba
      sed -i '1i#!/usr/bin/env python3' "$f"
      echo "Shebang añadido en $f"
    fi
  fi
done


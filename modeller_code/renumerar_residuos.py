#!/usr/bin/env python3
import sys


def renumerar_residuos(pdb_file):
    output_file = pdb_file.replace(".pdb", "_renum.pdb")

    with open(pdb_file, "r") as f_in, open(output_file, "w") as f_out:
        current_residue = None
        new_resnum = 0

        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                res_id = (line[17:20].strip(), line[21].strip(), line[22:26].strip())  # (resname, chain, old_resnum)

                if res_id != current_residue:
                    new_resnum += 1
                    current_residue = res_id

                # reemplazar número de residuo (col. 23-26)
                new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                f_out.write(new_line)
            else:
                f_out.write(line)

    print(f"Archivo renumerado creado: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python3 renumerar_residuos.py archivo.pdb")
        sys.exit(1)

    renumerar_residuos(sys.argv[1])


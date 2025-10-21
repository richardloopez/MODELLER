# MODELLER
Scripts to automatize deep search in modeller


PASOS PARA USAR EL CÓDIGO AUTOMÁTICO DE MODELLER

1.	Descargar el .pdb original 
2.	Recomendado: Eliminar átomos innecesarios (agua, iones...) y limpiar .pdb de cadenas innecesarias
3.	Guardar la estructura como .pdb en DS: importante para eliminar metadata
4.	Con Force: poner toda la estructura a modelar en la misma “Chain” (“A”)
5.	Con “renum-HETATM_residuos.py”: Si la estructura no incluye ADN u otros ligandos que se quieran incluir, lanzar el código como se indica dentro del mismo. Si lleva ligandos, se deben indicar al lanzar el código para que estos cambien sus cadenas de “ATOM” a “HETATM”. Nota: Esto hay que valorarlo porque las predicciones de estructura con ligandos a veces son malas, sobre todo si los ligandos están incompletos.
6.	Una vez tenemos el XXX_renum_HETATM.pdb, vamos a preparar el código de modeller.
7.	En config.py:	
a.	Seleccionar los parámetros deseados 
b.	Incluir la secuencias requeridas: Completa (UniProt) y del pdb. Para obtener la completa, esta se encuentra disponible en UniProt. Para obtener la del pdb, se puede buscar “pdb to fasta” en internet y salen conversores. Luego, comprobar que está correcta: debe tener el mismo número de aminoácidos que residuos tipo ATM en el XXX_renum_HETATM.pdb
c.	Si la estructura tiene ligandos HETATM, estos deben indicarse de la siguiente manera: 
8.	sequence_full = "QEKATPRNLPPRFQDGYYS/...........................--/-................------------"
9.	pdb_aa = "QEK/...........................--/-................------------"
En este caso había una cadena “A” con la proteína  y luego otras dos cadenas “C” y “B” con ADN. Se indican los residuos de ADN con “.” Las “-“ indican residuos faltantes en el pdb. No se tendrán en cuenta, pero mejor ponerlos por trazabilidad. Aun así, este caso no se recomienda. Sólo se deben utilizar HETATM si están completos o faltan las puntas. De lo contrario, pasarán cosas rarillas. 
d. A partir de la “sequence_full” (UniProt) se lanza una predicción de estructura secundaria para los loops. Web: https://bioinf.cs.ucl.ac.uk/psipred/ Vale con marcar la opción “PSIPRED 4.0 (Predict Secondary Structure) PSIPRED 4.0 (Predict Secondary Structure)” Se recomienda añadir correo. Los cálculos suelen tardar horas. Una vez terminado, se descarga el .zip y se guarda el .ss2.
e. Revisar las CPUs a usar: En “modeller_lanzador.sh” se detallan las CPUs que se pedirán. Luego, en “config.py” se configuran los trabajos que verdaderamente lanzará modeller.
	i. Dentro de un mismo nodo (no hay paralelización entre diferentes nodos), con todas las CPUs vacías, cada trabajo debería ocupar una CPU diferente. Sin embargo, esto no es óptimo. La máxima velocidad se consigue cuando todos los hilos de cada CPU están ocupados (1 trabajo ocupa 1 hilo), pero una CPU puede tener varios hilos (ej: 2). Entonces, si se ponen las mismas CPUs en “modeller_lanzador.sh” y en “config.py” se estarán usando a “la mitad de potencia” todos los nodos que se pidan. Por lo tanto, si cada CPU cuenta con 2 hilos, lo óptimo es pedir el doble de trabajos (.py) que de CPUs ) (.sh). Pero cuidado, ej: Si el nodo tiene 40 CPUs y nosotros hacemos (.sh = 20) (.py = 40) estaremos consumiendo recursos que SLURM “no sabe”. Como él piensa que sólo se están usando 20 CPUs, meterá trabajos de otras personas hasta alcanzar las 40 CPUs, por lo que estaremos incordiando a otros usuarios. Por lo tanto, la recomendación es pedir un nodo completo para nosotros y ahí usar el doble de hilos (trabajos) que CPUs. Ej:
Si tenemos: nodo01, CPUs = 96, threads/CPU = 2 la recomendación es: [ .sh = 96 | .py = 192 ] . Así tenemos el nodo completo para nosotros solos y dentro lo usamos como nos conviene.

f. Lanzar trabajo “sbatch modeller_lanzador.sh”
g. Importante: Casi al momento de empezar el trabajo se generan los archivos de alineamiento. Este suele funcionar muy bien, pero a veces comete errores. Es buena práctica revisarlo. Si está correcto, se deja continuar. Si no lo está, es recomendable para el script y relanzarlo con un alineamiento facilitado por el usuario (opción, en config.py: “USE_MANUAL_ALIGNMENT = False” se cambia a True). Abajo se facilitan los nombres. Nota: Se deben facilitar tanto el alineamiento normal como el cde. Los formatos y estructuras se pueden inducir a partir de los “fallidos”.
h. Nota: Se puede saber el progreso del script.
	i. Se generan todos los AUTO
	ii. Para cada AUTO que pase a LOOP, se generan todos los LOOPs. O sea, hasta que no se terminan los LOOPs del primero, no comienza el siguiente. 

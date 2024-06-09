# Proyecto Final:

## Análisis de Rendimiento de Dotplot

## Secuencial vs Paralelización

Este proyecto busca poner en práctica y evaluar la eficacia de tres métodos distintos para crear y analizar dotplots, una herramienta habitual en el ámbito de la bioinformática para contrastar secuencias de ADN o proteínas.

El propósito de este proyecto es no solo implementar tres variantes diferentes de dotplots, sino también llevar a cabo un análisis exhaustivo de su rendimiento. Los dotplots son una técnica esencial en bioinformática, permitiendo la comparación visual de secuencias biológicas para identificar similitudes y diferencias entre ellas. Esta investigación no solo se centra en la creación práctica de los dotplots, sino que también se adentra en la evaluación de su efectividad y utilidad en el contexto de análisis genómico y proteómico.

### Prerequisitos

El desarrollo del proyecto se llevó a cabo en Python 3.12.2, aprovechando la capacidad de computación en paralelo a través de las bibliotecas multiprocessing, mpi4py y pycuda. Para ejecutar el cálculo del dot-plot, se necesitan parámetros de entrada, incluidas las secuencias de referencia y consulta en formato fasta (o fna), los cuales deben ser especificados en la línea de comandos al momento de ejecutar el programa.

### Instalacion

Tener instalado python, para la posterior instalación de la slibrerias necesarias

A continuación, instale los paquetes de Python necesarios

```
pip install biopython
pip install numpy
pip install matplotlib
pip install mpi4py
pip install opencv-python
pip install tqdm==2.2.3

```

Para tener el proyecto de manera local, clone el repositorio en la terminar de su entorno de ejecución:

```
https://github.com/AlejandraMolinaAnduquia/DotplotADN_2024-1.git

```

### Ejecución

Para ejecutar el programa secuencial, ejecute el siguiente comando:

```
python proyecto.py --file1=./data/E_coli.fasta --file2=./data/Salmonella.fasta --sequential
```

Para ejecutar multiprocessing, ejecute el siguiente comando:

```
python proyecto.py --file1=./data/E_coli.fasta --file2=./data/Salmonella.fasta --multiprocessing
```

Para ejecutar mpi4py, ejecute el siguiente comando:

```
python proyecto.py --num_processes 1 2 3 4 --file1=./data/E_coli.fasta --file2=./data/Salmonella.fasta --mpi
```

Para ejecutar pycuda, ejecute el siguiente comando:

```

```
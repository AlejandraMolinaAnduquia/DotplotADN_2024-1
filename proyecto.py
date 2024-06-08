from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import argparse
import time
import cv2
from tqdm import tqdm

def read_fasta(file_name):
    sequences = []
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append(str(record.seq))
    return "".join(sequences)

def draw_dotplot(dotplot, fig_name='dotplot.svg'):
    plt.figure(figsize=(10, 10))
    plt.imshow(dotplot, cmap="Greys", aspect="auto")
    plt.xlabel("Secuencia 1")
    plt.ylabel("Secuencia 2")
    plt.savefig(fig_name)
    plt.show()

def dotplot_sequential(sequence1, sequence2):
    dotplot = np.empty((len(sequence1), len(sequence2)))
    for i in tqdm(range(len(sequence1))):
        for j in range(len(sequence2)):
            if sequence1[i] == sequence2[j]:
                if i == j:
                    dotplot[i, j] = 1
                else:
                    dotplot[i, j] = 0.7
            else:
                dotplot[i, j] = 0
    return dotplot

def apply_filter(matrix, path_image):
    kernel_diagonales = np.array([[1, -1, -1],
                                  [-1, 1, -1],
                                  [-1, -1, 1]])
    filtered_matrix = cv2.filter2D(matrix, -1, kernel_diagonales)
    normalized_matrix = cv2.normalize(filtered_matrix, None, 0, 127, cv2.NORM_MINMAX)
    threshold_value = 50
    _, thresholded_matrix = cv2.threshold(normalized_matrix, threshold_value, 255, cv2.THRESH_BINARY)
    cv2.imwrite(path_image, thresholded_matrix)
    cv2.imshow('Diagonales', thresholded_matrix)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--file1', dest='file1', type=str,
                        default=None, help='Query sequence in FASTA format')
    parser.add_argument('--file2', dest='file2', type=str,
                        default=None, help='Subject sequence in FASTA format')

    parser.add_argument('--sequential', action='store_true',
                        help='Ejecutar en modo secuencial')
    args = parser.parse_args()

    file_path_1 = args.file1
    file_path_2 = args.file2

    try:
        merged_sequence_1 = read_fasta(file_path_1)
        merged_sequence_2 = read_fasta(file_path_2)
    except FileNotFoundError as e:
        print("Archivo no encontrado, verifique la ruta")
        exit(1)

    Secuencia1 = merged_sequence_1[0:16000]
    Secuencia2 = merged_sequence_2[0:16000]

    if args.sequential:
        start_secuencial = time.time()
        dotplotSequential = dotplot_sequential(Secuencia1, Secuencia2)
        print(f"Tiempo de ejecución secuencial: {time.time() - start_secuencial}")
        draw_dotplot(dotplotSequential[:600, :600], fig_name="images/images_sequential/dotplot_secuencial.png")
        path_image = 'images/images_filter/dotplot_filter_sequential.png'  
        apply_filter(dotplotSequential[:600, :600], path_image)

if __name__ == "__main__":
    main()

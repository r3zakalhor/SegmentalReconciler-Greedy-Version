import os
import pandas as pd
import matplotlib.pyplot as plt
from PyPDF2 import PdfMerger
import numpy as np

# Replace 'your_folder_path' with the actual folder path containing your CSV files
folder_path = 'nb_dup_simphy_reconciliation_505'
folder_path2 = 'nb_dup_lca_reconciliation_505'
folder_path3 = 'nb_dup_greedy_reconciliation_505'
folder_path4 = 'nb_dup_ultragreedy_reconciliation_505'


# Get a list of CSV files in the folder
csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

# Create a PdfMerger object
pdf_merger = PdfMerger()

# Loop through each CSV file
for csv_file in csv_files:
    # Construct the full file path
    file_path = os.path.join(folder_path, csv_file)
    file_path2 = os.path.join(folder_path2, csv_file)
    file_path3 = os.path.join(folder_path3, csv_file)
    file_path4 = os.path.join(folder_path4, csv_file)

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(file_path)
    df2 = pd.read_csv(file_path2)
    df3 = pd.read_csv(file_path3)
    df4 = pd.read_csv(file_path4)

    # Assuming the CSV has two columns named 'species' and 'nb dup'
    x_values = df['species id']
    y_values = df['nb dup']
    x_values2 = df2['species id']
    y_values2 = df2['nb dup']
    x_values3 = df3['species id']
    y_values3 = df3['nb dup']
    x_values4 = df4['species id']
    y_values4 = df4['nb dup']
    # Create a line chart

    plt.plot(x_values, y_values, marker='', linestyle='-', label="simphy")
    plt.plot(x_values2, y_values2, marker='', linestyle='-', label="lca")
    plt.plot(x_values3, y_values3, marker='', linestyle='-', label="greedy")
    plt.plot(x_values4, y_values4, marker='', linestyle='-', label="ultragreedy")

    # Fit a polynomial of degree 2 (you can adjust the degree as needed)
    """
    coefficients = np.polyfit(x_values, y_values, 9)
    polynomial = np.poly1d(coefficients)
    x_poly = np.linspace(min(x_values), max(x_values), 100)
    y_poly = polynomial(x_poly)
    plt.plot(x_poly, y_poly, label="simphy")


    coefficients2 = np.polyfit(x_values2, y_values2, 9)
    polynomial2 = np.poly1d(coefficients2)
    x_poly2 = np.linspace(min(x_values2), max(x_values2), 100)
    y_poly2 = polynomial(x_poly2)
    plt.plot(x_poly2, y_poly2, label="lca")

    coefficients3 = np.polyfit(x_values3, y_values3, 9)
    polynomial3 = np.poly1d(coefficients3)
    x_poly3 = np.linspace(min(x_values3), max(x_values3), 100)
    y_poly3 = polynomial(x_poly3)
    plt.plot(x_poly3, y_poly3, label="greedy")

    coefficients4 = np.polyfit(x_values4, y_values4, 9)
    polynomial4 = np.poly1d(coefficients4)
    x_poly4 = np.linspace(min(x_values4), max(x_values4), 100)
    y_poly4 = polynomial(x_poly4)
    plt.plot(x_poly4, y_poly4, label="ultragreedy")
    """

    # Add labels and title
    plt.xlabel('Species id')
    plt.ylabel('NB Dup')
    plt.title(f'Line Chart from {csv_file}')
    # Add legend
    plt.legend()
    # Save the plot as a PDF file
    pdf_file_path = os.path.splitext(file_path)[0] + '.pdf'
    plt.savefig(pdf_file_path, format='pdf')

    # Close the current figure to start a new one for the next iteration
    plt.close()

    # Append the PDF file to the PdfMerger object
    pdf_merger.append(pdf_file_path)

# Output file path for the merged PDF
merged_pdf_path = os.path.join(folder_path, 'simphy_greedy_lca.pdf')

# Write the merged PDF to the output file
with open(merged_pdf_path, 'wb') as merged_output:
    pdf_merger.write(merged_output)

# Close the PdfMerger object
pdf_merger.close()

# Optional: Show a message when the merged PDF is created
print(f'Merged PDF saved at {merged_pdf_path}.')

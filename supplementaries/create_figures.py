import os
import pandas as pd
import matplotlib.pyplot as plt
from PyPDF2 import PdfMerger
import numpy as np

# Replace 'your_folder_path' with the actual folder path containing your CSV files
folder_path = 'nb_dup_species_sim4'

# Get a list of CSV files in the folder
csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

# Create a PdfMerger object
pdf_merger = PdfMerger()

# Loop through each CSV file
for csv_file in csv_files:
    # Construct the full file path
    file_path = os.path.join(folder_path, csv_file)

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(file_path)

    # Assuming the CSV has two columns named 'species' and 'nb dup'
    x_values = df['species id']
    y_values = df['nb dup']

    # Create a line chart
    plt.plot(x_values, y_values, marker='', linestyle='-')


    # Fit a polynomial of degree 2 (you can adjust the degree as needed)
    coefficients = np.polyfit(x_values, y_values, 9)
    polynomial = np.poly1d(coefficients)
    x_poly = np.linspace(min(x_values), max(x_values), 100)
    y_poly = polynomial(x_poly)
    plt.plot(x_poly, y_poly, label=f'Polynomial Trend Line (Degree {len(coefficients) - 1})')


    # Add labels and title
    plt.xlabel('Species id')
    plt.ylabel('NB Dup')
    plt.title(f'Line Chart from {csv_file}')

    # Save the plot as a PDF file
    pdf_file_path = os.path.splitext(file_path)[0] + '.pdf'
    plt.savefig(pdf_file_path, format='pdf')

    # Close the current figure to start a new one for the next iteration
    plt.close()

    # Append the PDF file to the PdfMerger object
    pdf_merger.append(pdf_file_path)

# Output file path for the merged PDF
merged_pdf_path = os.path.join(folder_path, 'merged_output.pdf')

# Write the merged PDF to the output file
with open(merged_pdf_path, 'wb') as merged_output:
    pdf_merger.write(merged_output)

# Close the PdfMerger object
pdf_merger.close()

# Optional: Show a message when the merged PDF is created
print(f'Merged PDF saved at {merged_pdf_path}.')

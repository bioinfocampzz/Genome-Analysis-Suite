# Genome Analysis Suite

A PyQt5-based desktop application for basic bioinformatics workflows. The app can load FASTA and FASTQ files, analyze DNA/RNA sequences, search for motifs, generate reverse complements and protein translations, and export results to CSV/Excel.

## Features

- Parse common sequence formats: FASTA, FNA, FASTQ, and FQ
- Validate IUPAC nucleotide sequences
- Calculate sequence length and GC content
- Generate reverse complements
- Translate sequences into amino acids
- Search for motifs, including overlapping matches
- Export analysis results to CSV and Excel
- Optional Plotly/WebEngine analytics dashboard

## Requirements

- Python 3.9 or later
- PyQt5
- PyQtWebEngine
- plotly
- openpyxl

Install the dependencies with:

```bash
pip install PyQt5 PyQtWebEngine plotly openpyxl
```

## Run the App

From the project folder, run:

```bash
python app.py
```

If your environment uses a different filename or launcher, use the script listed in this repository as the main entry point.

## Build a Windows Executable

The project includes a PyInstaller spec file. To build a packaged app:

```bash
pyinstaller app.spec
```

The build output will be created in the `dist/` folder.

## Project Structure

- `app.py` - main desktop application
- `app.spec` - PyInstaller build configuration
- `data/` - sample FASTA/FASTQ files for testing
- `genome_results.csv` - example output file
- `build/` - PyInstaller build artifacts

## Sample Data

Use the files in `data/` to test the application quickly:

- `data/example.fastq`
- `data/example.fna`
- `data/GCA_000195955.2_ASM19595v2_genomic.fna`

## Notes

- The app can run without Plotly WebEngine support, but the dashboard features will be limited if optional dependencies are missing.
- If you plan to upload this project to GitHub, consider adding a `.gitignore` file to exclude build artifacts such as `build/`, `dist/`, and Python cache files.

## License

Add your preferred license here before publishing the repository.

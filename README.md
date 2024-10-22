# waters2imzml

# Waters to imzML Converter

This script converts Waters raw data files to imzML format using ProteoWizard (through Docker). It takes a folder containing raw Waters datasets as input and converts each raw file to mzMLb format using ProteoWizard. Then, it extracts the x and y coordinates from the .inf file associated with each raw file and converts the mzMLb file to imzML format using pyteomics and pyimzml libraries.

## Requirements

- Python==3.11.8
- chardet==5.2.0
- pyimzML==1.5.3
- pyteomics==4.7.2
- ProteoWizard (installed as a Docker container)

## Usage

1. Install the required Python libraries:
    ```
    pip install -r requirements.txt
    ```

2. Install Docker Desktop from the official Docker website: https://www.docker.com/products/docker-desktop

3. Run the script:
    ```
    python waters_to_imzml.py
    ```

4. Follow the prompts to provide the path to the folder containing the raw Waters datasets, the polarity (positive or negative), and the position of the ID in the file name.

5. The script will convert each raw file to mzMLb format using ProteoWizard and then convert it to imzML format using pyteomics and pyimzml libraries. The converted imzML files will be saved in the 'imzml' folder within the raw data folder.


Shield: [![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg


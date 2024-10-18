from pyteomics import mzmlb
from pyimzml.ImzMLWriter import ImzMLWriter
import chardet
import pathlib
import os
import gc



class mzmlb_conv():
    """
    A class for converting mzMLb files to imzML format.

    Args:
        file_path (str): The path to the input mzMLb file.
        coord_tuple (tuple): A tuple containing the x and y dimensions of the image.
        output_file_path (str): The path to the output imzML file.
        polarity (str): The polarity of the spectra ('positive' or 'negative').

    Attributes:
        file_path (str): The path to the input mzMLb file.
        output_file_path (str): The path to the output imzML file.
        x (int): The x dimension of the image.
        y (int): The y dimension of the image.
        polarity (str): The polarity of the spectra ('positive' or 'negative').
        file (mzmlb.MzMLbFile): The mzMLb file object.
        coords (list): A list of coordinate tuples (x, y, z).

    Methods:
        _create_profile_spectra: Creates profile spectra from the mzMLb file.
        _convert_to_imzml: Converts the profile spectra to imzML format.
    """

    def __init__(self, file_path, coord_tuple, output_file_path, polarity):
        """
        Initializes the mzmlb_conv object.

        Args:
            file_path (str): The path to the input mzMLb file.
            coord_tuple (tuple): A tuple containing the x and y dimensions of the image.
            output_file_path (str): The path to the output imzML file.
            polarity (str): The polarity of the spectra ('positive' or 'negative').
        """

        self.file_path = file_path
        self.output_file_path = output_file_path
        self.x, self.y = coord_tuple
        self.polarity = polarity

        self.file = mzmlb.read(self.file_path)
        self.coords = []

        for i in range(self.x):
            for j in range(self.y):
                self.coords.append((i+1, j+1, 1))

        self._create_profile_spectra()
        self._convert_to_imzml()
        print('done.')
        gc.collect()

    def _create_profile_spectra(self):
        """
        Creates profile spectra from the mzMLb file using pyteomics.
        """

        self.spectra = []
        print('creating profile spectra...')
        for i in range(self.x*self.y):
            dict_pixel = self.file.get_by_index(i)
            self.spectra.append([dict_pixel['m/z array'], dict_pixel['intensity array'], self.coords[i]])

        gc.collect()

    def _convert_to_imzml(self):
        """
        Converts the profile spectra to imzML format using pyimzml.
        """

        print('converting to imzml...')
        with ImzMLWriter(self.output_file_path, polarity=self.polarity, spec_type='profile') as w:
            for mzs, intensities, coords in self.spectra:
                # writes data to the .ibd file
                w.addSpectrum(mzs, intensities, coords)
        gc.collect()

def get_encoding_type(file):
    """
    Get the encoding type of a file.

    Parameters:
    file (str): The path to the file.

    Returns:
    str: The encoding type of the file.
    """
    with open(file, 'rb') as f:
        result = chardet.detect(f.read())
    return result['encoding']


def get_coords(inf_file):
    """
    Retrieves the x and y coordinates from an input file.

    Args:
        inf_file (str): The path to the input file.

    Returns:
        tuple: A tuple containing the x and y coordinates.

    """
    encoding = get_encoding_type(inf_file)

    with open(inf_file, 'r', encoding=encoding) as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('DesiXLength'):
            x = float(line.split('\t')[5])
        if line.startswith('DesiXStep'):
            xstep = float(line.split('\t')[5])
        if line.startswith('DesiYLength'):
            y = float(line.split('\t')[5])
        if line.startswith('DesiYStep'):
            ystep = float(line.split('\t')[5])

    xcoord, ycoord = round(x / xstep), round(y / ystep)
    print('Image size:', xcoord, ' by ', ycoord)
    return (xcoord, ycoord)

raw_file_entries = []

def make_folders(raw_data_folder):
    """
    Creates necessary folders for processing raw data files.

    Parameters:
    raw_data_folder (str): The path to the raw data folder.

    Returns:
    tuple: A tuple containing the following elements:
        - nsplitfile (int): The number of splits in the raw data folder path.
        - mzmlb_folder (str): The path to the 'mzmlb' folder.
        - imzml_folder (str): The path to the 'imzml' folder.
        - pw_raw_data_folder (str): The path to the raw data folder with a prefix.
        - pw_output_folder (str): The path to the output folder with a prefix.
    """
    nsplitfile = len(raw_data_folder.split('/'))-1
  
    pw_raw_data_folder = raw_data_folder +':/' + raw_data_folder.split('/')[nsplitfile]
    mzmlb_folder = '/'.join([raw_data_folder, 'mzmlb/'])
    pw_output_folder = '/'.join(['',raw_data_folder.split('/')[nsplitfile], 'mzmlb'])

    try:
        os.makedirs(mzmlb_folder)
    except FileExistsError:
        pass

    imzml_folder = '/'.join([raw_data_folder, 'imzml'])

    try:
        os.makedirs(imzml_folder)
    except FileExistsError:
        pass

    return nsplitfile, mzmlb_folder, imzml_folder, pw_raw_data_folder, pw_output_folder


def run_conversion(raw_data_folder, nid, pmode):
    """
    Converts raw data files to imzML format using ProteoWizard.

    Args:
        raw_data_folder (str): The path to the folder containing the raw data files.
        nid (int): The index of the ID in the raw file name.
        pmode (str): The polarity mode for the conversion (e.g., 'positive', 'negative').

    Returns:
        None
    """

    nsplitfile, mzmlb_folder, imzml_folder, pw_raw_data_folder, pw_output_folder = make_folders(raw_data_folder)
    raw_file_entries = list(pathlib.Path(raw_data_folder).glob('*.raw'))

    for raw_file in raw_file_entries:
        print('Processing:', raw_file)
        raw_file = str(raw_file).split('/')[-1]
        pw_raw_file = '/'.join(['', raw_data_folder.split('/')[nsplitfile], raw_file])
        file_folder = '/'.join([raw_data_folder, raw_file])
        print('Make sure Docker is running.')
        print('Running ProteoWizard for conversion to mzMLb...')
        pw_command = 'docker run --rm -e WINEDEBUG=fixme-all+msgbox+relay -v ' + pw_raw_data_folder + ' chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ' + pw_raw_file + ' --mzMLb -o ' + pw_output_folder + ' >/dev/null 2>&1'
        # added ' >/dev/null 2>&1' at the end of command to suppress console output

        os.system(pw_command)  # how to add progress bar?

        id = raw_file.split('_')[nid]
        id_lookup = '*SAMPLE_' + id + '*.*'
        print(id_lookup)
        print(str(list(pathlib.Path(mzmlb_folder).glob(id_lookup))))

        mzmlb_file = str(list(pathlib.Path(mzmlb_folder).glob(id_lookup))[0])
        # try:
        #     mzmlb_file = str(list(pathlib.Path(mzmlb_folder).glob(id_lookup))[0])
        # except IndexError:
        #     print('Make sure Docker is working')
            

        inf_file = str(list(pathlib.Path(file_folder).glob('_*.inf'))[0])
        print('Accessing .inf file for coordinate extraction:', inf_file)  # or maybe just input it?/ or change to _extern.inf

        mzmlb_to_imzml = mzmlb_conv(mzmlb_file, get_coords(inf_file), imzml_folder + '/' + str(id) + '.imzml', polarity=pmode)
        gc.collect()

if __name__== '__main__':
    # convert to mzmlb using proteowizard
    # # check if spaces in names. replace with _

    raw_data_folder = input('Path to folder containing raw Waters datasets:')
    pmode = input('Input polarity positive/negative:')
    nid = int(input('Position of id in file name:'))
    run_conversion(raw_data_folder,nid,pmode)
   
    #  need to have docker desktop running
    docker_command = "docker run -d -p 8080:80 nginx"
    os.system(docker_command)

    # convert to mzmlb using proteowizard
    # # check if spaces in names. replace with _
   
    raw_data_folder = input('Path to folder containing raw Waters datasets:')
    pmode = input('Input polarity positive/negative:')
    nid = int(input('Position of id in file name:'))
    run_conversion(raw_data_folder,nid,pmode)
    
    
   
# close files
   
   

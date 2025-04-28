from pyteomics import mzmlb
from pyimzml.ImzMLWriter import ImzMLWriter
import numpy as np
import psutil
import os
import gc
import logging
from typing import List, Tuple
from joblib import Parallel, delayed, parallel_backend
import chardet
import pathlib

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MemoryMonitor:
    @staticmethod
    def check_memory(threshold_gb=14):
        """Check if available memory is below threshold"""
        available_gb = psutil.virtual_memory().available / (1024 ** 3)
        if available_gb < threshold_gb:
            logger.warning(f"Low memory warning: {available_gb:.1f}GB available")
            return False
        return True

    @staticmethod
    def force_garbage_collection():
        """Force garbage collection and memory cleanup"""
        gc.collect()
        if hasattr(psutil, 'Process'):
            try:
                process = psutil.Process(os.getpid())
                process.memory_info()
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass

class ChunkedSpectrumProcessor:
    def __init__(self, max_chunk_size_mb=500):
        self.max_chunk_size_mb = max_chunk_size_mb
        self.memory_monitor = MemoryMonitor()

    def process_chunk(self, start: int, end: int, file: object, coords: List[Tuple]) -> List:
        """Process a chunk of spectra with memory monitoring"""
        chunk_spectra = []
        batch_size = 5  # Process 10 spectra at a time
        
        try:
            for i in range(start, end, batch_size):
                batch_end = min(i + batch_size, end)
                
                # Check memory before processing batch
                if not self.memory_monitor.check_memory():
                    logger.warning("Low memory detected. Forcing garbage collection...")
                    self.memory_monitor.force_garbage_collection()
                
                batch_spectra = []
                for j in range(i, batch_end):
                    try:
                        dict_pixel = file.get_by_index(j)
                        # Use numpy's float32 to reduce memory usage
                        mz_array = dict_pixel["m/z array"].astype(np.float32) 
                        intensity_array = dict_pixel["intensity array"].astype(np.float32)
                        batch_spectra.append([mz_array, intensity_array, coords[j]])
                    except Exception as e:
                        logger.error(f"Error processing spectrum {j}: {str(e)}")
                        continue
                
                chunk_spectra.extend(batch_spectra)
                del batch_spectra
                gc.collect()
                
        except Exception as e:
            logger.error(f"Error in chunk {start}-{end}: {str(e)}")
            
        return chunk_spectra

class mzmlb_conv:
    def __init__(self, file_path, coord_tuple, output_file_path, polarity):
        self.file_path = file_path
        self.output_file_path = output_file_path
        self.x, self.y = coord_tuple
        self.polarity = polarity
        self.spectra = []
        self.processor = ChunkedSpectrumProcessor()
        
        # Calculate maximum chunk size based on available memory
        total_pixels = self.x * self.y
        self.chunk_size = min(1000, total_pixels // 4)  # Limit chunk size
        
        try:
            self.file = mzmlb.read(self.file_path)
            self.coords = [(i+1, j+1, 1) for i in range(self.x) for j in range(self.y)]
            self._create_profile_spectra()
            self._convert_to_imzml()
        finally:
            
            gc.collect()


    def _create_profile_spectra(self):
        """Creates profile spectra with memory-efficient parallel processing"""
        logger.info('Creating profile spectra...')
        
        try:
            # Use only 2 processes for large files to avoid memory issues
            n_jobs = 2
            total_pixels = self.x * self.y
            
            # Create smaller chunks
            chunk_indices = [
                (i, min(i + self.chunk_size, total_pixels))
                for i in range(0, total_pixels, self.chunk_size)]
            
            logger.info(f"Processing {len(chunk_indices)} chunks")
            
            # Process chunks in batches to control memory usage
            batch_size = 2
            for i in range(0, len(chunk_indices), batch_size):
                if not MemoryMonitor.check_memory():
                    logger.warning("Low memory detected. Waiting for cleanup...")
                    MemoryMonitor.force_garbage_collection()
                    
                current_chunks = chunk_indices[i:i + batch_size]
                
                with parallel_backend('loky', n_jobs=n_jobs):
                    results = Parallel(verbose=1)(
                        delayed(self.processor.process_chunk)(
                            start, end, self.file, self.coords
                        ) for start, end in current_chunks
                    )
                
                for chunk_result in results:
                    self.spectra.extend(chunk_result)
                
                del results
                gc.collect()
                
        except Exception as e:
            logger.error(f"Error in parallel processing: {str(e)}")
            raise

    def _convert_to_imzml(self):
        """Converts spectra to imzML format with memory monitoring"""
        logger.info('Converting to imzML...')
        
        try:
            with ImzMLWriter(self.output_file_path, polarity=self.polarity, spec_type='profile') as w:
                batch_size = 100  # Process spectra in smaller batches
                
                for i in range(0, len(self.spectra), batch_size):
                    if not MemoryMonitor.check_memory():
                        logger.warning("Low memory detected during conversion. Forcing cleanup...")
                        MemoryMonitor.force_garbage_collection()
                        
                    batch = self.spectra[i:i + batch_size]
                    for mzs, intensities, coords in batch:
                        w.addSpectrum(mzs, intensities, coords)
                    
                    # Clear processed batch
                    del batch
                    gc.collect()
                    
        except Exception as e:
            logger.error(f"Error during imzML conversion: {str(e)}")
            raise
        finally:
            # Clear all spectra data
            self.spectra.clear()
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


def run_conversion(raw_data_folder, nid, pmode, mzmlb_status):
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
        if mzmlb_status:
           
            print('Make sure Docker is running.')
            print('Running ProteoWizard for conversion to mzMLb...')
            pw_command = 'docker run --rm -e WINEDEBUG=fixme-all+msgbox+relay -v ' + pw_raw_data_folder + ' chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ' + pw_raw_file + ' --mzMLb -o ' + pw_output_folder + ' >/dev/null 2>&1'
            # added ' >/dev/null 2>&1' at the end of command to suppress console output
            # print(pw_command)
            os.system(pw_command)  # how to add progress bar?

        id = raw_file.split('_')[nid]
        id_lookup = '*_' + id + '_*.*'
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
    import sys
  
    
    raw_data_folder = str(sys.argv[1])
    pmode = str(sys.argv[2])
    nid = int(sys.argv[3])
    mzmlb_status = eval(sys.argv[4])

    print(mzmlb_status)

    # raw_data_folder = input('Path to folder containing raw Waters datasets:')
    # pmode = input('Input polarity positive/negative:')
    # nid = int(input('Position of id in file name:'))
    run_conversion(raw_data_folder,nid,pmode, mzmlb_status)
   
    #  need to have docker desktop running
    docker_command = "docker run -d -p 8080:80 nginx"
    os.system(docker_command)
  
    # close files
   
   

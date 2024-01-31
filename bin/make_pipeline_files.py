import argparse
import glob
import sys
import os
from datetime import date

class DirectoryFiles:
    def __init__(self, top_dir, pip_dir):
        """Initialize class

        Inputs -
            top_dir: str - top directory to build off
            pip_dir: str - pipeline files go here
        """
        self.top_dir = os.path.abspath(top_dir) # top directory to build from
        self.pip_dir = os.path.abspath(pip_dir) # pipeline files

        self.log_dir = {'node': os.path.abspath(os.path.join(top_dir, 'log'))}
        self.ana_dir = {'node': os.path.abspath(os.path.join(top_dir, 'analysis'))}
        self.ben_dir = {'node': os.path.abspath(os.path.join(top_dir, 'benchmark'))}

    def get_analysis_node_dir(self):
        return self.ana_dir['node']

    def get_benchmark_node_dir(self):
        return self.ben_dir['node']

    def get_log_node_dir(self):
        return self.log_dir['node']

    def get_analysis_lane_dir(self, lane):
        return self.ana_dir[lane]

    def get_benchmark_lane_dir(self, lane):
        return self.ben_dir[lane]

    def get_log_lane_dir(self, lane):
        return self.log_dir[lane]

    def create_base_directories(self):
        # Create top level directory
        try:
            os.mkdir(self.top_dir)
        except FileExistsError:
            check = input(f'WARNING: {self.top_dir} already exists! Overwrite (Y) or Skip (n)? [n]/Y ')
            if check == 'Y':
                print(f'WARNING: You chose to overwrite the pipeline files!')
            else:
                print(f'MESSAGE: Skipping {self.top_dir}')
                return False

        # Create pipeline files directory
        os.makedirs(self.pip_dir, exist_ok=True)

        # Create processed data directory and subdirectories
        os.makedirs(self.get_log_node_dir(), exist_ok=True)
        os.makedirs(self.get_analysis_node_dir(), exist_ok=True)
        os.makedirs(self.get_benchmark_node_dir(), exist_ok=True)

        return True

    def create_lane_directories(self, lane):
        # Add lane directories to dicts
        self.log_dir[lane] = os.path.abspath(os.path.join(self.get_log_node_dir(), lane))
        self.ana_dir[lane] = os.path.abspath(os.path.join(self.get_analysis_node_dir(), lane))
        self.ben_dir[lane] = os.path.abspath(os.path.join(self.get_benchmark_node_dir(), lane))

        # Create directories
        os.makedirs(self.get_log_lane_dir(lane), exist_ok=True)
        os.makedirs(self.get_analysis_lane_dir(lane), exist_ok=True)
        os.makedirs(self.get_benchmark_lane_dir(lane), exist_ok=True)

# Helper function to format date
def format_date():
    return date.today().strftime('%d %B %Y')   

# Helper function to print message and exit
def print_error_and_exit(error):
    print(error)
    sys.exit(1)

# Helper function to read file contents into one string
def get_file_contents(fname):
    with open(fname, 'r') as fh:
        template = fh.read()
    return template

# Helper function to get template config file
def get_config_template():
    return get_file_contents('../config/config.yaml')

# Helper function to get template slurm file
def get_submit_template():
    return get_file_contents('submit_workflow_template.slurm')

def parse_args():
    """Parse command line arguments.

    Inputs -
        None
    Returns -
        argparse.ArgumentParser.parse_args()
    """
    parser = argparse.ArgumentParser(
        description = 'Create pipeline input files from a file with a list of names and directories'
    )

    parser.add_argument(
        'directory_file',
        help='file with name prefixes and directories'
    )

    return parser.parse_args()

def parse_input_file(fname):
    """Parse input file with prefixes and directories

    Inputs -
        fname: str - filename to parse
    Returns -
        list: elements of (prefix, directory)
    """
    out = []
    with open(fname, 'r') as fh:
        try:
            for l in fh.readlines():
                line = l.strip().split('\t')

                prefix = line[0]
                direct = line[1]

                # Catch early if input directory doesn't exist
                if not os.path.exists(direct):
                    raise FileNotFoundError(f'ERROR: {direct} does not exist. Please check and rerun.')

                out.append((prefix, direct))
        except FileNotFoundError as e:
            print_error_and_exit(e)

    return out

def setup_base_directory_structure(prefix):
    """Create directory structure for new run

    Inputs -
        prefix: str - name of top directory to create for new processing run
    Returns -
        DirectoryFiles if created else None
    """
    top_dir = f'../results/{prefix}'
    pip_dir = f'{top_dir}/pipeline_files'

    dir_files = DirectoryFiles(top_dir, pip_dir)
    created = dir_files.create_base_directories()
    if not created:
        return None

    return dir_files

def format_well_name(well):
    """Format name of a well to work as input to platetools

    Inputs -
        well: str - well name to format (should be a single letter followed by either a 1 or 2 digit number)
    Returns -
        str: formatted well name
    """
    if len(well) == 3:
        return well
    elif len(well) == 2:
        letter = well[0]
        number = int(well[1])
        return f'{letter}{number:02d}'
    else:
        print_error_and_exit(f'ERROR: malformed well name input: {well}')

def get_id_and_lane(fname):
    """Parse FASTQ name to extract sample ID and lane information

    Inputs -
        fname: str - FASTQ file name
    Return -
        tuple: (str, str) - (sample ID, lane)
    """
    pieces = [x for x in fname.split('_') if 'plate' not in x]

    # Get index of lane element, this will be the position before the read name
    # Handling it this way avoids relying on the sample ID always being "Well_CellId"
    lane_index = None
    if 'R1' in pieces:
        lane_index = pieces.index('R1') - 1
    elif 'R2' in pieces:
        lane_index = pieces.index('R2') - 1
    else:
        print_error_and_exit('ERROR: cannot retrieve sample ID and lane info from {fname}. Missing R1 and R2')

    # Correctly format well
    pieces[0] = format_well_name(pieces[0])

    return ('_'.join(pieces[0:lane_index]), pieces[lane_index])

def create_samplesheet(dname):
    """Write samplesheet(s) from FASTQ files in a directory

    Inputs -
        dname: str - directory name to search
    Returns -
        list: samplesheet file names created
    """
    # Get read 1 FASTQs (read 2 name based on read 1 name)
    fastqs = glob.glob(f'{dname}/*_R1_001.fastq.gz')

    # If no files match the pattern, break early
    if len(fastqs) == 0:
        print_error_and_exit(f'ERROR: {dname} has no files that match the pattern *_R1_001.fastq.gz')

    # Multiple lanes may exist in one directory, split those up here
    lanes = {}
    try:
        for fpath in fastqs:
            path, f1 = os.path.split(fpath)
            f2 = f1.replace('_R1_001.fastq.gz', '_R2_001.fastq.gz')

            # Check read 2 exists
            if not os.path.exists(os.path.join(path, f2)):
                raise FileNotFoundError(f'ERROR: malformed read 2 file name - {f2} could not be found.')

            sample_id, lane = get_id_and_lane(f1)

            try:
                lanes[lane].append((sample_id, f1, f2))
            except KeyError:
                lanes[lane] = [(sample_id, f1, f2)]
    except FileNotFoundError as e:
        print_error_and_exit(e)

    # Write out samplesheets and collect file names for returning from function
    samplesheets = []
    for lane, files in lanes.items():
        sample_fname = f'{lane}_samples.tsv'
        samplesheets.append(sample_fname)

        with open(sample_fname, 'w') as fh:
            fh.write('sample\tfq1\tfq2\n')
            for id, fq1, fq2 in files:
                fh.write(f'{id}\t{fq1}\t{fq2}\n')

    return samplesheets

def create_config_file(template, dname, snames, dir_files):
    """Create config file from template in config/ directory

    Inputs -
        template: str - config file template
        dname: str - directory where FASTQs live
        snames: str - samplesheet name
        dir_files: DirectoryFiles - paths where output should be written
    Returns -
        list: config file names created
    """
    # Do some modifications ahead that don't rely on samplesheet
    template = template.replace(
        'STORMqc Configuration Template',
        'AUTO-GENERATED CONFIG! DO NOT MODIFY!\n#\n# STORMqc Pipeline'
    )
    template = template.replace('CREATE_DATE', format_date())
    template = template.replace('FASTQ_DIR', dname)

    configs = []
    for sname in snames:
        lane = sname.replace('_samples.tsv', '')
        config_fname = f'{lane}_config.yaml'
        configs.append(config_fname)

        current = template
        current = current.replace('ANALYSIS_DIR', dir_files.get_analysis_lane_dir(lane))
        current = current.replace('BENCHMARK_DIR', dir_files.get_benchmark_lane_dir(lane))
        current = current.replace('LOG_DIR', dir_files.get_log_lane_dir(lane))
        current = current.replace('SAMPLESHEET', f'{os.getcwd()}/{sname}')

        with open(config_fname, 'w') as fh:
            fh.write(current)

    return configs

def create_submit_file(template, cnames, dir_files):
    """Create slurm submit file from template in bin/ directory

    Inputs -
        template: str - submit file template
        lname: str - log directory path
        cnames: str - config file names
        dir_files: DirectoryFiles - paths where output should be written
    Returns -
        list: submit file names created
    """
    # Do some modifications ahead that don't rely on samplesheet
    template = template.replace(
        'STORMqc SLURM Submit Template',
        'AUTO-GENERATED SUBMIT SCRIPT! DO NOT MODIFY!\n#\n# STORMqc Submit Script'
    )
    template = template.replace('CREATE_DATE', format_date())

    submit_files = []
    for cname in cnames:
        lane = cname.replace('_config.yaml', '')
        submit_fname = f'{lane}_submit.slurm'
        submit_files.append(submit_fname)

        current = template
        current = current.replace('__CONFIG_YAML__', f'{dir_files.pip_dir}/{cname}')
        current = current.replace('__LOG_DIRECTORY__', dir_files.get_log_lane_dir(lane))

        with open(submit_fname, 'w') as fh:
            fh.write(current)

    return submit_files

if __name__ == '__main__':
    # Command line arguments
    args = parse_args()

    # Get config file template before we start moving directories
    CONFIG_TEMPLATE = get_config_template()

    # Get slurm submit file template before we start moving directories
    SUBMIT_TEMPLATE = get_submit_template()

    # Set a place to come back to before creating a new set of directories
    HOMEBASE = os.getcwd()

    # Parse input file to get directories to create pipeline files for
    dirs = parse_input_file(args.directory_file)

    # Loop over directories creating all necessary directories and files
    for p, d in dirs:
        # Create directories and switch to new pipeline files directory
        # Try creating directories and switch to new pipeline if desired
        # Can't create lane directories since we haven't processed the input directory yet
        dir_files = setup_base_directory_structure(p)
        if dir_files is None:
            continue
        os.chdir(dir_files.pip_dir)

        # Make samplesheets
        samplesheets = create_samplesheet(d)

        # Now that we've processed the input directory, we can create lane subdirectories
        lanes = [s.replace('_samples.tsv', '') for s in samplesheets]
        for lane in lanes:
            dir_files.create_lane_directories(lane)

        # Make config files
        configs = create_config_file(CONFIG_TEMPLATE, d, samplesheets, dir_files)

        # Make submit scripts
        submits = create_submit_file(SUBMIT_TEMPLATE, configs, dir_files)

        os.chdir(HOMEBASE)

        print(f'Sucessfully created directory structure and pipelines files for {p}')

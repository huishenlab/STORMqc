import argparse
import glob
import gzip
import sys
import os
from datetime import date

# Retrieved from StackOverflow:
#     https://stackoverflow.com/questions/287871/how-do-i-print-colored-text-to-the-terminal
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class DirectoryFiles:
    def __init__(self, top_dir, pip_dir, overwrite):
        """Initialize class

        Inputs -
            top_dir: str - top directory to build off
            pip_dir: str - pipeline files go here
            overwrite: str - how to handle overwriting of existing directories
        """
        self.overwrite = overwrite
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
            if self.overwrite == 'ask':
                print_warning(f'{self.top_dir} already exists!')
                check = input('Overwrite (Y) or Skip (n)? [n]/Y ')
                if check == 'Y':
                    print_warning('You chose to overwrite the pipeline files!')
                else:
                    print_message(f'Skipping {self.top_dir}')
                    return False
            elif self.overwrite == 'always':
                    print_warning('You chose to overwrite the pipeline files!')
            else:
                print_message(f'Skipping {self.top_dir}')
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

# Helper function to print error message and exit
def print_error_and_exit(error):
    print(f'{bcolors.FAIL}[ERROR] {error}{bcolors.ENDC}')
    sys.exit(1)

# Helper function to print warning
def print_warning(warning):
    print(f'{bcolors.WARNING}[WARNING] {warning}{bcolors.ENDC}')

# Helper function to print message
def print_message(message):
    print(f'{bcolors.OKGREEN}[MESSAGE] {message}{bcolors.ENDC}')

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

# Helper function to get template slurm file for rerunning make_plots
def get_rerun_template():
    return get_file_contents('rerun_quality_control_plots_template.slurm')

# Helper function to check minimum read count met
def has_min_read_count(fname, min):
    # Assume gzipped files that are 10kb or larger are close enough to the 100 limit
    # TODO: This may need to be changed from a hardcoded value if we find most people
    #       up the limit over 100 reads
    if os.stat(fname).st_size > 10000:
        return True

    with gzip.open(fname, 'rb') as fh:
        n_lines = sum(1 for line in fh)

    if n_lines < 4*min:
        return False

    return True

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
        '-m', '--min-read-count',
        default=100,
        type=int,
        help='minimum number of reads required in FASTQ file'
    )

    parser.add_argument(
        '-O', '--overwrite',
        choices = ['ask', 'always', 'never'],
        default='ask',
        help='how to handle existing experiment directories (ask [default]: ask for confirmation | always: always overwrite | never: never overwrite)'
    )

    parser.add_argument(
        '-V', '--visualize-config',
        default='default_visualization_limits.yaml',
        help='YAML file controlling visualization limits'
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
                    raise FileNotFoundError(f'{direct} does not exist. Please check and rerun.')

                out.append((prefix, direct))
        except FileNotFoundError as e:
            print_error_and_exit(e)
        except IndexError:
            print_error_and_exit(f'malformed input file: {fname}')

    return out

def setup_base_directory_structure(prefix, overwrite):
    """Create directory structure for new run

    Inputs -
        prefix: str - name of top directory to create for new processing run
        overwrite: str - how to handle overwriting of existing directories
    Returns -
        DirectoryFiles if created else None
    """
    top_dir = f'../results/{prefix}'
    pip_dir = f'{top_dir}/pipeline_files'

    dir_files = DirectoryFiles(top_dir, pip_dir, overwrite)
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
        print_error_and_exit(f'malformed well name input: {well}')

def split_name_and_remove_unwanted_tags(fname):
    """Split filename while also removing unwanted tags.

    Inputs -
        fname: str - file name
    Returns -
        list, split on '_'
    """
    pieces = fname.split('_')

    problem_tags = ['plate', 'Plate']
    for problem in problem_tags:
        pieces = [x for x in pieces if problem not in x]

    return pieces

def get_id_and_lane(fname):
    """Parse FASTQ name to extract sample ID and lane information

    Inputs -
        fname: str - FASTQ file name
    Return -
        tuple: (str, str) - (sample ID, lane)
    """
    pieces = split_name_and_remove_unwanted_tags(fname)

    # Get index of lane element, this will be the position before the read name
    # Handling it this way avoids relying on the sample ID always being "Well_CellId"
    lane_index = None
    if 'R1' in pieces:
        lane_index = pieces.index('R1') - 1
    elif 'R2' in pieces:
        lane_index = pieces.index('R2') - 1
    else:
        print_error_and_exit(f'cannot retrieve sample ID and lane info from {fname}. Missing R1 and R2')

    # Correctly format well
    pieces[0] = format_well_name(pieces[0])

    return ('_'.join(pieces[0:lane_index]), pieces[lane_index])

def create_samplesheet(dname, min_reads):
    """Write samplesheet(s) from FASTQ files in a directory

    Inputs -
        dname: str - directory name to search
        min_reads: int - minimum number of reads to include file
    Returns -
        list: samplesheet file names created
    """
    # Get read 1 FASTQs (read 2 name based on read 1 name)
    fastqs = glob.glob(f'{dname}/*_R1_001.fastq.gz')

    # If no files match the pattern, break early
    if len(fastqs) == 0:
        print_error_and_exit(f'{dname} has no files that match the pattern *_R1_001.fastq.gz')

    # Multiple lanes may exist in one directory, split those up here
    lanes = {}
    too_low = {}
    try:
        for fpath in fastqs:
            path, f1 = os.path.split(fpath)
            f2 = f1.replace('_R1_001.fastq.gz', '_R2_001.fastq.gz')

            # Check read 2 exists
            if not os.path.exists(os.path.join(path, f2)):
                raise FileNotFoundError(f'malformed read 2 file name - {f2} could not be found.')

            sample_id, lane = get_id_and_lane(f1)

            # Check for minimum read count
            # Read 1 and Read 2 should have the same number of reads, so only check read 1
            has_min = has_min_read_count(fpath, min_reads)
            if not has_min:
                try:
                    too_low[lane].append(sample_id)
                except KeyError:
                    too_low[lane] = [sample_id]
            else:
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

    # Write out list of samples with too few reads
    for lane, files in too_low.items():
        output_fname = f'{lane}_too_few_reads.txt'

        # Only write output file if samples exist with too few reads
        if len(files) > 0:
            with open(output_fname, 'w') as fh:
                fh.write('sample\n')
                for id in files:
                    fh.write(f'{id}\n')

    return samplesheets

def create_config_file(template, dname, snames, dir_files, resource_paths, viz_config):
    """Create config file from template in config/ directory

    Inputs -
        template: str - config file template
        dname: str - directory where FASTQs live
        snames: str - samplesheet name
        dir_files: DirectoryFiles - paths where output should be written
        resource_paths: dict - paths to various resources
        viz_config: str - path to visualization config file
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
    template = template.replace('STAR_INDEX', resource_paths['star_idx'])
    template = template.replace('RRNA_PATH', resource_paths['rrna_path'])
    template = template.replace('MITO_PATH', resource_paths['mito_path'])
    template = template.replace('ERCC_PATH', resource_paths['ercc_path'])
    template = template.replace('ANNOT_SPACE', resource_paths['annt_path'])
    template = template.replace('VIZ_LIMITS', viz_config)

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

def create_submit_file(template, file_tag, cnames, dir_files):
    """Create slurm submit file from template in bin/ directory

    Inputs -
        template: str - submit file template
        file_tag: str - tag to include between lane name and file extension [(lane)_(tag).slurm]
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
        submit_fname = f'{lane}_{file_tag}.slurm'
        submit_files.append(submit_fname)

        current = template
        current = current.replace('__CONFIG_YAML__', f'{dir_files.pip_dir}/{cname}')
        current = current.replace('__LOG_DIRECTORY__', dir_files.get_log_lane_dir(lane))

        with open(submit_fname, 'w') as fh:
            fh.write(current)

    return submit_files

def check_resource_files_exist(base_dir):
    """Checks if needed input files have been created already.

    Inputs -
        None
    Returns -
        bool, exits if missing
    """
    SOLUTION = 'Trying running resources/gather_rna_resources.slurm'

    if not os.path.exists(f'{base_dir}/resources/star_index_human'):
        print_error_and_exit(f'STAR index for the human genome is missing. {SOLUTION}')
    if not os.path.exists(f'{base_dir}/resources/star_index_mouse'):
        print_error_and_exit(f'STAR index for the mouse genome is missing. {SOLUTION}')

    if not os.path.exists(f'{base_dir}/resources/GRCh38_ERCC_enhancers_rmsk.merged.sorted.bed.gz'):
        print_error_and_exit(f'Annotated space BED file is missing. {SOLUTION}')

    if not os.path.exists(f'{base_dir}/resources/gene_ids.hg38.rrna.txt'):
        print_error_and_exit(f'human rRNA gene IDs missing. {SOLUTION}')
    if not os.path.exists(f'{base_dir}/resources/gene_ids.hg38.mito.txt'):
        print_error_and_exit(f'human mtDNA gene IDs missing. {SOLUTION}')
    if not os.path.exists(f'{base_dir}/resources/gene_ids.hg38.ercc.txt'):
        print_error_and_exit(f'human ERCC gene IDs missing. {SOLUTION}')

    if not os.path.exists(f'{base_dir}/resources/gene_ids.mm10.rrna.txt'):
        print_error_and_exit(f'mouse rRNA gene IDs missing. {SOLUTION}')
    if not os.path.exists(f'{base_dir}/resources/gene_ids.mm10.mito.txt'):
        print_error_and_exit(f'mouse mtDNA gene IDs missing. {SOLUTION}')
    if not os.path.exists(f'{base_dir}/resources/gene_ids.mm10.ercc.txt'):
        print_error_and_exit(f'mouse ERCC gene IDs missing. {SOLUTION}')

    return True

if __name__ == '__main__':
    # Command line arguments
    args = parse_args()
    if args.overwrite == 'ask':
        print_message('You have chosen to confirm overwriting of all directories that already exist!')
    elif args.overwrite == 'always':
        print_message('You have chosen to overwrite all directories that already exist!')
    else:
        print_message('You have chosen to ignore all directories that already exist!')

    # Get config file template before we start moving directories
    CONFIG_TEMPLATE = get_config_template()

    # Get slurm submit file template before we start moving directories
    SUBMIT_TEMPLATE = get_submit_template()
    RERUN_TEMPLATE = get_rerun_template()

    # Set a place to come back to before creating a new set of directories
    HOMEBASE = os.getcwd()

    # Check if resources files have already been created
    # TODO: Allow user to select either human or mouse
    check_resource_files_exist(f'{HOMEBASE}/..')
    resource_paths = {
        'star_idx': os.path.abspath(f'{HOMEBASE}/../resources/star_index_human'),
        'rrna_path': os.path.abspath(f'{HOMEBASE}/../resources/gene_ids.hg38.rrna.txt'),
        'mito_path': os.path.abspath(f'{HOMEBASE}/../resources/gene_ids.hg38.mito.txt'),
        'ercc_path': os.path.abspath(f'{HOMEBASE}/../resources/gene_ids.hg38.ercc.txt'),
        'annt_path': os.path.abspath(f'{HOMEBASE}/../resources/GRCh38_ERCC_enhancers_rmsk.merged.sorted.bed.gz'),
    }

    # Check visualization config exists
    if not os.path.exists(args.visualize_config):
        print_error_and_exit(f'visualization config file does not exist: {args.visualize_config}')
    viz_config = os.path.abspath(args.visualize_config)

    # Parse input file to get directories to create pipeline files for
    dirs = parse_input_file(args.directory_file)

    # Loop over directories creating all necessary directories and files
    for p, d in dirs:
        # Create directories and switch to new pipeline files directory
        # Try creating directories and switch to new pipeline if desired
        # Can't create lane directories since we haven't processed the input directory yet
        dir_files = setup_base_directory_structure(p, args.overwrite)
        if dir_files is None:
            continue
        os.chdir(dir_files.pip_dir)

        # Make samplesheets
        samplesheets = create_samplesheet(d, args.min_read_count)

        # Now that we've processed the input directory, we can create lane subdirectories
        lanes = [s.replace('_samples.tsv', '') for s in samplesheets]
        for lane in lanes:
            dir_files.create_lane_directories(lane)

        # Make config files
        configs = create_config_file(CONFIG_TEMPLATE, d, samplesheets, dir_files, resource_paths, viz_config)

        # Make submit scripts
        submits = create_submit_file(SUBMIT_TEMPLATE, 'submit', configs, dir_files)
        reruns = create_submit_file(RERUN_TEMPLATE, 'rerun', configs, dir_files)

        os.chdir(HOMEBASE)

        print_message(f'Sucessfully created directory structure and pipelines files for {p}')

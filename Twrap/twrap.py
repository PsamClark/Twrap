#!/usr/bin/env python

import datetime
import multiprocessing as mp
import os
import re
import time
import shutil
import subprocess
import sys
from subprocess import PIPE, Popen
import pandas as pd
import glob
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

class InputError(Exception):
    """Error raised for arguments"""

    pass


class FileError(Exception):
    """Error raised for files"""
    pass


class SeqError(Exception):
    """Error raised for for sequences"""
    pass


class colour:
    """Colours for progress bar"""

    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def time_stamp():
    """
    Returns the date and time.
    """
    return '%s' % ( datetime.datetime.fromtimestamp(
        time.time()).strftime('%Y-%m-%d %H:%M:%S'))


def run_jobs(jobs_to_do, cpus=2):
    """
    Runs jobs that need to be submitted for processing.

    Args:
        cpus (int): number of CPUs to use
    """
    running_jobs = {}
    total_jobs = len(jobs_to_do)

    # print initial progress
    print(time_stamp(), 'Total jobs to do: %d' %
          (total_jobs))
    print(time_stamp(), '%d%% of jobs completed.' % ( int(
        (total_jobs - len(jobs_to_do)) / float(total_jobs) * 100)))

    # run through the jobs to do
    while len(jobs_to_do) > 0:
        finished_jobs = []

        # identify finished jobs
        for job in running_jobs:
            alive = job.is_alive()
            if not alive:
                finished_jobs.append(job)

        # remove finished jobs from running jobs list
        if len(finished_jobs) > 0:
            print(time_stamp(), '%d%% of jobs completed' %
                  (int((total_jobs - len(jobs_to_do)) /
                       float(total_jobs) * 100)))

        # remove complete jobs from job list
        for job in finished_jobs:
            running_jobs.pop(job)

        # try start more jobs
        jobs_to_run = cpus - len(running_jobs)

        if jobs_to_run > 0:
            for job in range(jobs_to_run):
                if len(jobs_to_do) > 0:
                    new_job = jobs_to_do.pop(0)
                    new_job.start()
                    running_jobs[new_job] = 0
        else:
            time.sleep(1)

    while len(running_jobs) > 0:
        finished_jobs = []
        for job in running_jobs:
            alive = job.is_alive()
            if not alive:
                finished_jobs.append(job)
        for job in finished_jobs:
            running_jobs.pop(job)

    # print progress
    print(time_stamp(), '%d%% of jobs completed' %
          (int((total_jobs - len(jobs_to_do)) / float(total_jobs) * 100)))


class Twrap(object):
    """Motif detection class

    Attributes:
    _seq_file (str): FASTA file containing sequences.
    seqs (:obj: 'list' of :obj: 'str'): Sequences under investigation.
    seq_names (:obj: 'list' of :obj: 'str'): Names of sequences under
    investigation.
    _shp_file (str): File containing shape data.
    _shp (:obj: 'list' of :obj: 'float'): Shape data profile.
    _sl_method (str): Stemloop/loop determining method.
    _win_size (int): Size of folding window for scanning.
    _win_inc (int): Increment window scanned moved by in each iteration.
    _ecut (float): Stemloop energy cutoff filter.
    _fcut (int): Stemloop/loop frequency cutoff filter.
    _mb_file (str): File containing m&b values.
    _mb_vals (np.ndarray): Array containing m&b values.
    _nuc_type (str): Type of nucleotides in sequence.
    _param_file (str): Parameter file employed in folding.
    _temp_seq (str): Temporary sequence storage.
    _temp_name (str): Temporary sequence name storage.
    _shp_id (int): Identifier for m&b values used.
    _temp_shp (:obj: 'list' of :obj: 'float'): Temporary shape data storage.
    _sl_data ('list'): Stemloop/loop information storage.
    """

    def __init__(self, seqs_in):

        # sequence attributes
        self._seq_file = ''
        self.seqs = []
        self.seq_names = []
        self.seq_input(seqs_in)

        # shape related attributes
        self._shp_file = None
        self._shp = []
        self._mb_file = None
        self._mb_vals = []

        # folding attributes
        self._sl_method = 'stemloops'
        self._win_size = None
        self._win_inc = None
        self._ecut = None
        self._fcut = 0
        self._nuc_type = 'RNA'
        self._param_file = ''

        # temporary storage
        self._temp_seq = ''
        self._temp_name = ''
        self._shp_id = ''
        self._temp_shp = []
        self._sl_data = []

    def _input_checker(self):
        """Check the arguments input are accceptable"""

        # check folding window size
        if self._win_size is not None:
            if not isinstance(self._win_size, int):
                raise InputError('Window size must be integer!')
            elif not np.isfinite(self._win_size):
                raise InputError('Window size must be finite!')
            elif self._win_size < 0:
                raise InputError('Window size must be positive!')

        # check window sliding increment
        if self._win_inc is not None:
            if not isinstance(self._win_inc, int):
                raise InputError('Window increment must be integer!')
            elif not np.isfinite(self._win_inc):
                raise InputError('Window increment must be finite!')
            elif self._win_inc < 0:
                raise InputError('Window increment must be positive!')

        # specify stemloop/loop methods possible
        methods = ['stemloops', 'loops']
        if not isinstance(self._sl_method, str):
            raise InputError('sl method must be string!')

        elif self._sl_method not in methods:
            raise InputError('sl method not found!')

        # check energy cutoff
        if self._ecut is not None:
            if not isinstance(self._ecut, float):
                raise InputError('energy cutoff must be float!')
            elif not np.isfinite(self._ecut):
                raise InputError('energy cutoff must be finite!')

        # check frequency cutoff
        if self._fcut is not None:
            if not isinstance(self._fcut, int):
                raise InputError('Frequency cutoff must be integer!')
            elif not np.isfinite(self._win_inc):
                raise InputError('Frequency cutoff must be finite!')
            elif self._win_inc < 0:
                raise InputError('Frequency cutoff must be positive!')

        # check shape file
        if self._shp_file is not None:
            if not isinstance(self._shp_file, str):
                raise InputError('Shape file must be string or None!')

            elif not os.path.isfile(sys.path[0]+'/' + self._shp_file):
                raise FileError(
                    'Shape file "' +
                    self._shp_file +
                    '" not found!')

        # specify nucleotide types that can be folded
        nuc_types = ['RNA', 'DNA']

        # nucleotide type checks
        if not isinstance(self._nuc_type, str):
            raise InputError('Nucleotide type must be string!')

        elif self._nuc_type.upper() not in nuc_types:
            raise InputError('Nucleotide type not found!')

        # check m&b file
        if self._mb_file is not None:
            if not isinstance(self._mb_file, str):
                raise InputError('m&b file must be string or None!')

            elif not os.path.isfile(sys.path[0]+'/' + self._mb_file):
                raise FileError('m&b file "' + self._mb_file + '" not found!')

    def _file_dir_pop(self):

        # delete and recreate fold_output dir
        if os.path.exists(sys.path[0]+'/fold_output'):
            shutil.rmtree(sys.path[0]+'/fold_output')
        os.makedirs(sys.path[0]+'/fold_output')

        # create stemloop/loop output dir
        try:
            os.mkdir(sys.path[0]+'/'+self._sl_method + '_output')
        except OSError:
            pass

        # create shape data output file
        try:
            os.mkdir(sys.path[0]+'/shp_output')
        except OSError:
            pass

        # print sequence output dir
        try:
            os.mkdir(sys.path[0]+'/seq_output')
        except OSError:
            pass

        # if m&b file exists populat value array
        if self._mb_file is not None:
            mb_vals = pd.read_csv(sys.path[0]+'/'+self._mb_file)
            self._mb_vals = mb_vals.values.astype(float)

        # if shape file exists populate data
        if self._shp_file is not None:

            self._shp = [float(i.strip('\n')) for i in 
                open(sys.path[0]+'/'+self._shp_file, 'r').readlines()]

        # setup parameter file
        param_files = {
            'RNA': 'rna_turner1999.dat',
            'DNA': 'dna_santalucia.dat'}

        self._param_file = param_files[self._nuc_type]

    def fold(self, win_size=None, win_inc=None, sl_method='stemloops',
             ecut=None, fcut=0, shp_file=None, nuc_type='RNA', mb_file=None):
        """Fold sequence wrapper.

        Args:
            win_size (int): Size of folding window for scanning.
            win_inc (int): Increment window scanned moved by in each iteration.
            sl_method (str): Stemloop/loop determining method.
            ecut (float): Stemloop energy cutoff filter.
            fcut (int): Stemloop/loop frequency cutoff filter.
            shp_file (str): File containing shape data.
            nuc_type (str): Type of nucleotides in sequence.
            mb_file (str): File containing m&b values.
        """

        # assign arguments to attributes
        self._win_size = win_size
        self._win_inc = win_inc
        self._sl_method = sl_method
        self._ecut = ecut
        self._fcut = fcut
        self._shp_file = shp_file
        self._nuc_type = nuc_type
        self._mb_file = mb_file

        # check attributes
        self._input_checker()

        # create files and directories
        self._file_dir_pop()

        # populate processes
        if self._shp_file is None:
            processes = [mp.Process(target=self.fold_threaded, args=(i,))
                         for i in range(len(self.seqs))]
        else:
            processes = []
            for j in range(len(self.seqs)):
                processes = [mp.Process(target=self.fold_threaded, args=(
                    j, i,)) for i in range(len(self._mb_vals))]
                processes += sub_processes

        # run processes
        run_jobs(processes, cpus=10)

    def fold_threaded(self, i, mb=None):
        """Multithread processing of folding.

        Args:
            i (int): Specifies the sequence to be folded.
            mb (int): Specifies the m&b values to use.
        """

        # if m&bs used follow this path
        if mb is not None:
            self.fold_wrapper(i, mb)
            self.sl_reader_mb(i, mb)

        # else follow this path
        else:
            self.fold_wrapper(i)
            self.sl_reader(i)

            for filename in Path(sys.path[0]+'/shp_output').glob(self.seq_names[i]+'*'):
                filename.unlink()

            for filename in Path(sys.path[0]+'/seq_output').glob(self.seq_names[i]+'*'):
                filename.unlink()

    def seqfile_restructure(self):
        """Restructure sequence files so that they can be read by Tfold"""

        # if sequence file doesnt already exist continue
        if not os.path.exists(sys.path[0]+'/seq_output/' + self._temp_name):

            # create and write info into seq_file
            seqf_out = open(sys.path[0]+'/seq_output/' + self._temp_name, 'w')
            seqf_out.write('>' + self._temp_name + '\n')
            i = 0

            # write sequence 60 characters per line
            while i < len(self._temp_seq):
                seqf_out.write(self._temp_seq[i:i + 60] + '\n')
                i += 60
            seqf_out.close()

    def seq_input(self, seqin=''):
        """ Loads the sequence to be searched given a file path
            to a fasta file containing the sequence.
        Args:
            seq_in (str): Sequence input.
        """

        # if seqin is file
        if isinstance(seqin, str):
            if not os.path.isfile(sys.path[0]+'/' + seqin):
                raise FileError('Sequence file "' + seqin + '" not found!')

            # iterate through sequence file and extract sequence information
            with open(sys.path[0]+'/'+seqin, 'r') as handle:
                seq_recs = SeqIO.parse(handle, "fasta")

                # check file in FASTA format
                if not any(seq_recs):
                    raise FileError('File not in FASTA format!')

            # get sequence records
            seq_recs = SeqIO.parse(seqin, 'fasta')

        # if seqin is list of SeqRecord objects
        elif isinstance(seqin, list) and \
                all(isinstance(s, SeqRecord) for s in seqin):
            seq_recs = seqin

        # if seqin neither raise error
        else:
            raise SeqError('Incorrect sequence input!')

        # read through sequence records
        for record in seq_recs:
            seq = str(record.seq).upper()

            # check sequence records actually have sequences
            if len(seq) == 0:
                raise SeqError('Empty sequence detected in record: ' +
                               record.id)

            # check exclusivity of T and U use
            elif seq.count('T') > 0 and seq.count('U') > 0:
                raise SeqError('Sequence ' + record.id +
                               ' contains both U and T bases!')

            self.seq_names.append(str(record.id).replace(' ', ''))
            self.seqs.append(str(record.seq).upper())

        # if no seqs present raise error
        if len(self.seqs) == 0:
            raise FileError('File given had no sequences in it!')

    def shape_writer(self):
        "write shape profiles to file for use my Tfold"""

        # check file doesn't already exist
        if not os.path.exists(sys.path[0]+'/shp_output/' + \
            self._shp_id + '.shp'):
            f = open(sys.path[0]+'/shp_output/' + \
                self._shp_id + '.shp', 'w')

            # write data to file
            for i in self._temp_shp[:-1]:
                f.write(str(i) + '\n')
            f.write(str(self._temp_shp[-1]))

            f.close()

    def fold_wrapper(self, i, mb=None):
        """Folding wrapper algorithm
        Args:
            i (int): Indicator of sequence under investigation
            mb (int): Indicator of m&b values used
        """

        # get sequence info
        seq_name = self.seq_names[i]
        seq_in = self.seqs[i]

        # if global folding requesting
        if self._win_size is None:

            # populating temp seq info
            self._temp_name = seq_name
            self._temp_seq = seq_in

            # restructure sequence in file for use by Tfold
            self.seqfile_restructure()

            # create temporary shape data array
            if self._shp_file is not None:
                self._temp_shp = self._shp

            # if no shape given make pseudo shape array
            else:
                self._temp_shp = [-999] * len(seq_in)

            # write shape to file for use by Tfold
            self.shape_writer()

            # run TFold
            self.run_tfold(mb)

        # if local folding requested
        else:

            # iterate across genome with sliding window
            for j in range(-self._win_size + self._win_inc, len(seq_in) +
                           self._win_size, self._win_inc):

                # if start of window negative just take the start as zero
                if j < 0:
                    jb = 0
                    jt = j + self._win_size

                else:
                    jb = j
                    jt = j + self._win_size

                # create temp sequence for window
                self._temp_seq = seq_in[jb:jt]
                if self._shp_file is not None:
                    self._temp_shp = self._shp[jb:jt]
                else:
                    self._temp_shp = [-999] * len(self._temp_seq)

                # write shape to file for Tfold
                if len(self._temp_seq) > 4:
                    self._shp_id = seq_name+'_'+str(j) + '_' + str(jt)
                    self.shape_writer()

                    self._temp_name = seq_name + '_' + str(j) + '_' + str(jt)

                    # write seq to file for Tfold
                    self.seqfile_restructure()

                    # run Tfold
                    self.run_tfold(mb)

    def run_tfold(self, mb=None):
        """Run Tfold

        Args:
            mb (int): Indicator of m&b values used.
        """

        # if no shape data used
        if self._shp_file is None:

            # create bash command
            bashCommand = '%sTwrap/execs/Tfold.x -P %s \
                 -p %d -i %s -o %s -s %d -shp %s' % (sys.path[-1],
                    sys.path[-1]+'Twrap/params/'+self._param_file,
                    1000, sys.path[0]+'/seq_output/'+self._temp_name,
                    sys.path[0]+'/fold_output/' + self._temp_name + '.folds', 
                    123456789,
                    sys.path[0] + '/shp_output/' + self._shp_id + '.shp')

            # create and communicate process
            process = subprocess.Popen(
                bashCommand.split(), stdout=subprocess.PIPE)
            process.communicate()

        else:

            # create bash command
            bashCommand = '%s/Twrap/execs/Tfold.x -P %s \
                 -p %d -i %s -o %s -s %d -shp %s -m %f -b %f' % \
                (sys.path[-1], sys.path[-1]+'/Twrap/params/'+self._param_file,
                 1000, sys.path[0]+'/seq_output/'+self._temp_name,
                 sys.path[0]+'/fold_output/' + self._temp_name + '_' +
                 str(self._mb_vals[mb, 0]) + '_' +
                 str(self._mb_vals[mb, 1]) + '.folds', 123456789,
                 sys.path[0]+'/shp_output/' + self._shp_id + '.shp',
                 self._mb_vals[mb, 0], self._mb_vals[mb, 1])

            # create and communicate process
            process = subprocess.Popen(
                bashCommand.split(), stdout=subprocess.PIPE)
            process.communicate()

    def sl_reader(self, i):
        """Stemloop/loop reader.

        Args:
            i (int): Indicator of sequence used.
        """

        # get fold files
        fold_files = glob.glob(sys.path[0]+'/fold_output/' + self.seq_names[i] + '*')

        # open stemloop/loop output file
        data_file = open(
            sys.path[0]+'/' +
            self._sl_method +
            '_output/' +
            self.seq_names[i] +
            '_' +
            self._sl_method +
            '.csv',
            'w')

        # if stemloops to be calculated
        if 'stem' in self._sl_method:

            # write file headers
            data_file.write(
                'sl_pos,sl_seq,sl_fold,loop_pos,loop_seq,freq,energy\n')
        else:
            # write headers for loop only search
            data_file.write('loop_pos,loop_seq,freq\n')

        # create stemloop/loop data storage
        self._sl_data.append([{} for i in range(len(self.seqs[i]))])

        # run through fold files
        for ffilen in fold_files:
            ffile = open(ffilen, 'r')
            folds = ffile.read().split('\n')
            seq = folds[0]

            # get fold file ids
            fold_ids = ffilen.split('.folds')[0].split('_')

            # specify start position for fold
            if self._win_size is None:
                start_pos = 0
            else:
                start_pos = int(fold_ids[-2])

            # extract stemloops/loops from folds
            self.sl_populator(folds[1:], start_pos, seq, ffilen)

        # read through sl/l data and filter
        for iii in range(len(self._sl_data[-1])):
            for key, value in self._sl_data[-1][iii].items():
                sl_ids = key.split('_')

                # filter by frequency
                if value >= self._fcut:

                    if 'stem' in self._sl_method:

                        # calculate sl energy
                        cmd = sys.path[-1]+'/Twrap/execs/Eval.x -s ' + \
                            sl_ids[-2] + ' -f "' + sl_ids[-1] + '" -P ' + \
                            sys.path[-1]+'/Twrap/params/'+self._param_file
                        e_log = subprocess.check_output(
                            cmd, shell=True).decode()
                        if self._ecut is not None:

                            # filter by energy if cutoff given
                            if float(e_log) <= self._ecut:

                                # write to file
                                data_file.write(
                                    '%d,%s,%s,%d,%s,%d,%f\n' %
                                    (int(
                                        sl_ids[2]), sl_ids[3], sl_ids[4],
                                        int(sl_ids[0]) + start_pos,
                                        sl_ids[1], value, float(e_log)))
                        else:
                            data_file.write(
                                '%d,%s,%s,%d,%s,%d,%f\n' %
                                (int(
                                    sl_ids[2]), sl_ids[3], sl_ids[4],
                                    int(sl_ids[0]) + start_pos, sl_ids[1],
                                    value, float(e_log)))

                    else:
                        data_file.write('%d,%s,%d\n' %
                                        (int(sl_ids[0]), sl_ids[1], value))

        data_file.close()

    def sl_reader_mb(self, i, mb):
        """ Stemloop/loop reader when m&b's employed.

        Args:
            i (int): Indicator of sequence under investigation
            mb (int): Indicator of m&b values used
        """

        # get fold files
        fold_files0 = glob.glob(sys.path[0]+'/fold_output/' + self.seq_names[i] + '*')

        # filter down fold files based on m&b values
        fold_files = [cc for cc in fold_files0 if str(
            self._mb_vals[mb, 0]) + '_' + str(self._mb_vals[mb, 1]) in cc]

        # create sl/l data output file
        data_file = open(sys.path[0]+'/' +
                         self._sl_method +
                         '_output/' +
                         self.seq_names[i] +
                         '_' +
                         str(self._mb_vals[mb, 0]) +
                         '_' +
                         str(self._mb_vals[mb, 1]) +
                         '_' +
                         self._sl_method +
                         '.csv', 'w')
        if 'stem' in self._sl_method:
            # create header for sl
            data_file.write(
                'sl_pos,sl_seq,sl_fold,loop_pos,loop_seq,freq,energy\n')
        else:
            # create header for l
            data_file.write('loop_pos,loop_seq,freq\n')

        # create sl/l storage
        self._sl_data.append([{} for i in range(len(self.seqs[i]))])

        # read through fold files
        for ffilen in fold_files:
            ffile = open(ffilen, 'r')
            folds = ffile.read().split('\n')
            seq = folds[0]

            # get fold ids
            fold_ids = ffilen.split('.folds')[0].split('_')

            # get start position of fold
            if self._win_size is None:
                start_pos = 0
            else:
                start_pos = int(fold_ids[-4])

            # get sl/l data from folds
            self.sl_populator(folds[1:], start_pos, seq, ffilen)

        # read through sl/l data
        for ii in range(len(self._sl_data[-1])):
            for key, value in self._sl_data[-1][ii].items():
                sl_ids = key.split('_')

                # filter by frequency
                if value >= self._fcut:

                    if 'stem' in self._sl_method:

                        # get energy of stemloop
                        cmd = sys.path[-1]+'/Twrap/execs/Eval.x -s ' + \
                            sl_ids[-2] + ' -f "' + sl_ids[-1] + '" -P ' + \
                            sys.path[-1]+'/Twrap/params/'+self._param_file

                        e_log = subprocess.check_output(
                            cmd, shell=True).decode()
                        if self._ecut is not None:

                            # filter by energy cutoff
                            if float(e_log) <= self._ecut:

                                # write data to file
                                data_file.write(
                                    '%d,%s,%s,%d,%s,%d,%f\n' %
                                    (int(sl_ids[2]), sl_ids[3], sl_ids[4],
                                        int(sl_ids[0]), sl_ids[1], value,
                                        float(e_log)))
                        else:
                            data_file.write(
                                '%d,%s,%s,%d,%s,%d,%f\n' %
                                (int(sl_ids[2]), sl_ids[3], sl_ids[4],
                                    int(sl_ids[0]), sl_ids[1], value,
                                    float(e_log)))

                    else:
                        data_file.write('%d,%s,%d\n' %
                                        (int(sl_ids[0]), sl_ids[1], value))

        data_file.close()

    def sl_populator(self, folds, start, seq, ffilen):
        """Stemloop/loop data extractor

        Args:
            folds (:obj: 'list' of :obj: 'str'): Folds in vienna notation.
            start (int): Start position of folds in genome.
            seq (str): Sequence thats been folded.
            ffilen (str): Fold file name.
        """

        # iterate through folds
        for fold in folds:

            # find loops in fold
            for structure_match in re.finditer(
                    '[\\(\\.]+\\(\\.+?\\)[\\)\\.]+', fold):
                start_structure = structure_match.start()
                fold_structure = structure_match.group()

                # get open and close bracket positions
                opens = [
                    i for i in range(
                        len(fold_structure)) if fold_structure[i] == '(']
                closes = [
                    i for i in range(
                        len(fold_structure)) if fold_structure[i] == ')']
                # clip if more open or closes discovered
                if len(opens) < len(closes):

                    clip_ind = closes[len(opens) - 1]

                    fold_structure = fold_structure[:clip_ind + 1]
                elif len(opens) > len(closes):

                    clip_ind = opens[-len(closes)]

                    fold_structure = fold_structure[clip_ind:]

                    start_structure = start_structure + clip_ind

                # get new open and close positions
                opens = [
                    i for i in range(
                        len(fold_structure)) if fold_structure[i] == '(']
                closes = [
                    i for i in range(
                        len(fold_structure)) if fold_structure[i] == ')']

                # clip sl fold
                fold_structure = fold_structure[opens[0]:closes[-1] + 1]

                start_structure += opens[0]
                new_start = 0

                # get clipped sl seq
                subseq = seq[start_structure:start_structure +
                             len(fold_structure)]

                # get new start position
                if start > 0:
                    new_start = start

                # recalculate opens and closes
                opens = [
                    i for i in range(
                        len(fold_structure)) if fold_structure[i] == '(']
                closes = [
                    i for i in range(
                        len(fold_structure)) if fold_structure[i] == ')']
                # calculate sl and l start
                sl_start = start_structure + new_start + 1
                l_start = sl_start + opens[-1]

                # get loop sequence
                loopseq = subseq[opens[-1] + 1:closes[0]]

                # create or append data storage unit corresponding to sl/l
                if 'stem' in self._sl_method:
                    sl_id = str(l_start) + '_' + loopseq + '_' + \
                        str(sl_start) + '_' + subseq + '_' + fold_structure
                else:
                    sl_id = str(l_start) + '_' + loopseq

                if sl_id in list(self._sl_data[-1][l_start].keys()):

                    self._sl_data[-1][l_start][sl_id] += 1
                else:
                    self._sl_data[-1][l_start][sl_id] = 1

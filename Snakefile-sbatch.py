#!/usr/bin/env python
import sys
import subprocess
import os
import math
import ipdb
import errno
from snakemake.utils import read_job_properties

def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


class SnakeJob:
    """Snakemake can generate bash scripts that can be sumbitted by a
    scheduler.  This class reads the bash script and stores the number of the
    rule, name of bash file and the supplied input files."""
    def __init__(self, snakebashfile, dependencies=None):
        self.scriptname = snakebashfile
        job_properties = read_job_properties(snakebashfile)
        self.rule = job_properties['rule']
        self.ifiles = job_properties['input']
        self.ofiles = job_properties['output']
        if dependencies == None or len(dependencies) < 1:
            self.dependencies = None
        else:
            # expects snakemake like list of numbers
            self.dependencies = dependencies
            assert len(self.dependencies) >= 1

class UndefinedJobRule(Exception):
    """Exception in case an sbatch job has no defined resource usage in the
    code."""
    def __init__(self, msg):
        self.msg = msg


class SnakeJobSbatch(SnakeJob):
    # Change this to the path of the sbatch_job wrapper script
    sbatch_job_path = "./sbatch_job"
    proj_name = "b2014206"

    def __init__(self, snakebashfile, dependencies=None):
        SnakeJob.__init__(self, snakebashfile, dependencies)
        if self.dependencies == None:
            self.dep_str = ''
        else:
            self.dep_str = '-d ' + ','.join(["afterok:%s" % d for d in self.dependencies])

    def schedule(self):
        """Schedules a snakemake job with sbatch and determines resource usage
        based on input files."""
        # create the output directory, so slurm output can go there
        #make_dir(os.path.dirname(os.path.abspath(self.ofiles[0])))

        run_locally = False
        print(self.rule, file=sys.stderr)
        if self.rule == 'htseq_tophat_cutadapt':
            # Dummy rule, does not need any time
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '05',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}

            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {hours}:{minutes}:00 \
                    -J {job_name} {sbatch_job_path} \
                    '{script_name}'""".format(**attributes)
        elif self.rule == 'htseq_tophat_transcriptome_only_cutadapt':
            # Dummy rule, does not need any time
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '05',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}

            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {hours}:{minutes}:00 \
                    -J {job_name} {sbatch_job_path} \
                    '{script_name}'""".format(**attributes)
        elif self.rule == 'rna_seqc_tophat_cutadapt':
            # Dummy rule, does not need any time
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '05',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}

            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {hours}:{minutes}:00 \
                    -J {job_name} {sbatch_job_path} \
                    '{script_name}'""".format(**attributes)
        elif self.rule == 'qc':
            # c.a. 30M-40M reads c.a. 100bp long, longest job about 6 minutes
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '30',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'cutadapt':
            # c.a. 30M-40M reads c.a. 100bp long, longest job about 15 minutes
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '01',
                    'minutes': '00',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'tophat_index':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '01',
                    'minutes': '00',
                    'p': 'node',
                    'N': '1',
                    'n': '16',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'tophat':
            # c.a. 30M-40M reads c.a. 100bp long, longest job about 6 hours
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '07',
                    'minutes': '00',
                    'p': 'node',
                    'N': '1',
                    'n': '16',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'tophat_to':
            # c.a. 30M-40M reads c.a. 100bp long, longest job about 6 hours
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '07',
                    'minutes': '00',
                    'p': 'node',
                    'N': '1',
                    'n': '16',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'bam_index':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '20',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'bam_XS_index':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '20',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'bam_no_XS_index':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '20',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'count':
            # c.a. 30M-40M reads c.a. 100bp long, longest job about 35 min
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '01',
                    'minutes': '00',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'merge_count':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '00',
                    'minutes': '05',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'read_group':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '01',
                    'minutes': '00',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'reorder':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '01',
                    'minutes': '00',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)
        elif self.rule == 'rnaseqc':
            attributes = {
                    'dep_str': self.dep_str,
                    'days': '0',
                    'hours': '01',
                    'minutes': '00',
                    'p': 'core',
                    'N': '1',
                    'n': '1',
                    'job_name': "snakemake_{0}".format(self.rule),
                    'sbatch_job_path': self.sbatch_job_path,
                    'script_name': self.scriptname,
                    'proj_name': self.proj_name}
            sbatch_cmd = """sbatch {dep_str} -A {proj_name} -p {p} -N {N} -n {n} -t {days}-{hours}:{minutes}:00 \
                            -J {job_name} {sbatch_job_path} \
                            '{script_name}'""".format(**attributes)                                                        
        else:
            raise UndefinedJobRule('Undefined resource usage %s' % (self.rule))
            return 2

        print(sbatch_cmd, file=sys.stderr)
        popenrv = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).communicate()

        if not run_locally:
            try:
                print("%i" % int(popenrv[0].split()[-1]))
            except ValueError:
                print("Not a submitted job: %s" % popenrv[0])
                sys.exit(2)


if __name__ == '__main__':
    sj = SnakeJobSbatch(sys.argv[-1], sys.argv[1:-1])
    print(sys.argv[-1], sys.argv[1:-1], file =sys.stderr)
    try:
        sj.schedule()
    except UndefinedJobRule as err:
        print(err.msg, file=sys.stderr)
        sys.exit(2)

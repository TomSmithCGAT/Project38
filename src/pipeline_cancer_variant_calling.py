################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
Pipeline cancer_variant_calling
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline template.

Overview
========

The purpose of this pipeline is to take BAM files from cancer exome
sequencing and identify variants using Mutect.

Matched control/tumour samples are expected. Identifiers for the matched 
samples are given in the pipeline.ini file, 
e.g Control=Control, Tumour=Adeno
"ID-Control.1.bwa.bam"
"ID-Adeno.1.bwa.bam"

The pipeline will then run Mutect as specified in pipeline.ini






Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============




The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_cancer_variant_calling.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_cancer_variant_calling.tgz
   tar -xvzf pipeline_cancer_variant_calling.tgz
   cd pipeline_cancer_variant_calling
   python <srcdir>/pipeline_cancer_variant_calling.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# define some tracks if needed
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob("*.ini" ), "(\S+).ini" )


###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################

BAMS = glob.glob("*.bam")
reference_genome = PARAMS["mutect_reference"]


@follows(mkdir("realign.dir"))
@transform(BAMS,
           regex("(\S+).bam"),
           add_inputs(reference_genome),
           r"realign.dir/\1.grp")
def runBaseRecalibrator(infiles, outfile):

    infile, reference = infiles
    dbsnp = PARAMS["mutect_dbsnp"]

# need to unload java before runnning GATK as it now runs on java version 7

    statement = '''module unload apps/java/jre1.6.0_26;
    java -Xmx4g -jar
    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
    -T BaseRecalibrator
    -R %(reference)s
    -knownSites %(dbsnp)s
    -I %(infile)s
    -o %(outfile)s
    ''' % locals()

    P.run()


@transform(BAMS,
           regex("(\S+).bam"),
           add_inputs(reference_genome),
           r"realign.dir/\1_recalibrated.bam")
def recalibrateBAM(infiles, outfile):

    infile, reference = infiles

# need to unload java before runnning GATK as it now runs on java version 7

    statement = '''module unload apps/java/jre1.6.0_26;
    java -Xmx4g -jar
    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
   -T PrintReads
   -R %(reference)s
   -I %(infile)s
   -BQSR %(report)s
   -o %(outfile)s
    ''' % locals()

    P.run()


@transform(recalibrateBAM,
           suffix("(\S+).recalibrated.bam"),
           add_inputs(reference_genome),
           ".intervals")
def createRealignIntervals(infiles, outfile):

    infile, reference = infiles

# need to unload java before runnning GATK as it now runs on java version 7

    statement = '''module unload apps/java/jre1.6.0_26;
    java -Xmx4g -jar
    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
    -T RealignerTargetCreator
    -R %(reference)s
    -I %(infile)s
    -o %(outfile)s
    ''' % locals()

    P.run()


@follows(createRealignIntervals)
@transform(createRealignIntervals,
           suffix(".intervals"),
           add_inputs(reference_genome),
           r".realigned.bam")
def runRealign(infiles, outfile):

    intervals, reference = infiles
    infile = os.path.basename(intervals).replace(".intervals", ".bam")

# need to unload java before runnning GATK as it now runs on java version 7

    statement = '''module unload apps/java/jre1.6.0_26;
    java -Xmx4g -jar
    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
    -T IndelRealigner
    -R %(reference)s
    -I %(infile)s
    -targetIntervals %(intervals)s
    -o %(outfile)s''' % locals()

    P.run()


@follows(mkdir("mutect.dir/"))
@transform(runRealign,
           regex("realign.dir/(\S+)-Control(\S+).realigned.bam"),
           add_inputs(reference_genome),
           r"mutect.dir/\1.vcf")
def runMutect(infiles, outfile):

    infile, reference = infiles

    infile_tumour = infile.replace(
        "Control", PARAMS["mutect_tumour"])

    basename = P.snip(outfile, ".vcf")
    call_stats_out = basename + ".call_stats.out"
    coverage_wig_out = basename + ".coverage.wig"
    mutect_log = basename + ".mutect.log"

    cosmic, dbsnp, = (
        PARAMS["mutect_cosmic"],
        PARAMS["mutect_dbsnp"])

    statement = '''java -Xmx2g -jar
    /ifs/apps/bio/muTect-1.1.4/muTect-1.1.4.jar
    --analysis_type MuTect
    --reference_sequence %(reference)s
    --cosmic %(cosmic)s
    --dbsnp %(dbsnp)s
    --input_file:normal %(infile)s
    --input_file:tumor %(infile_tumour)s
    --out %(call_stats_out)s
    --coverage_file %(coverage_wig_out)s
    --vcf %(outfile)s
    --fraction_contamination 0.3
    --min_qscore 20
    --gap_events_threshold 1
    > %(mutect_log)s''' % locals()

    P.run()


@transform(runMutect,
           suffix(".vcf"),
           ".annotated.vcf")
def annotateVCF(infile, outfile):

    statement = '''java -Xmx4G -jar /ifs/apps/bio/snpEff-3.3-dev/snpEff.jar
                 -c /ifs/apps/bio/snpEff-3.3-dev/snpEff.config
                 hg19 %(infile)s > %(outfile)s''' % locals()

    P.run()

###################################################################
###################################################################
###################################################################
# primary targets
###################################################################


@follows(annotateVCF)
def full():
    pass


@follows(runBaseRecalibrator)
def partial():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish():
    '''publish report and data.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit(P.main(sys.argv))

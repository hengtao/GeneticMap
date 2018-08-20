#!/usr/bin/env python3

import re,time,os,sys
import traceback
import argparse
import json
import subprocess
import configparser
from functools import reduce
import logging
from datetime import datetime

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def Makedirs(dir_list):
    for everydir in dir_list:
        if not os.path.exists(everydir):
            os.makedirs(everydir)

def Traceback(traceback_recordfile):
    with open(traceback_recordfile, 'w+') as f:
        traceback.print_exc(file = f)
        f.flush()

def LogRecord(logger, loginfo):
    now_time = datetime.now()
    logger.info(loginfo)
    return now_time, logger

def Check_software(logger, software_path):
    if os.path.exists(software_path):
        logger.debug("Choose software:" + software_path + "!")
    else:
        output = os.popen('which ' + software_path).read()
        if output:
            software_temp = output.split("\n")[0]
            if os.path.exists(software_temp):
                software_path = software_temp
                logger.debug("Choose software:" + software_path + "!")
        else:
            location = os.popen("locate " + software_path).read()
            print(location)
            if location:
                software_path = location.split("\n")[0]
                print(software_path)
                if not os.path.exists(software_path):
                    logger.error("Can't locate the " + software_path + "!")
                    exit(1)
    return software_path

def Cmd_and_Time(cmd, logger):
    starttime, logger = LogRecord(logger, cmd)
    p = subprocess.Popen(cmd, shell = True)
    p.wait()
    if(p.returncode != 0):
        logger.exception("COMMAND fail:----\n{cmd}".format(cmd = cmd))
        sys.exit(1)
    else:
        spend_time = datetime.now() - starttime
        days       = spend_time.days
        totalseconds = spend_time.seconds
        hours      = int(int(totalseconds)/(60*60))
        minites    = int((int(totalseconds)%3600)/60)
        seconds    = int(totalseconds)%60
        logger.info("Total time for this step is: {days} days, {hours} hours, {minites} minites, {seconds} seconds".format(days = days, hours = hours, minites = minites, seconds = seconds))
    return spend_time

def VersionConvert():
	
def main(args):
    logfile  = args.logfile

    ## Import logger, 获取logger实例
    logger    = logging.getLogger(__name__)
    ## 指定logger输出格式
    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    #logger.setLevel(level = logging.INFO)
    handler = logging.FileHandler(logfile)
    handler.formatter = formatter
    ## 控制台日志
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.formatter = formatter
    ## 为logger添加日志处理器
    logger.addHandler(handler)
    logger.addHandler(console_handler)
    ## 指定日志的最低输出级别，默认为WARN
    logger.setLevel(level = logging.INFO)

    gatk     = Check_software(logger, "GenomeAnalysisTK.jar")
    vcftools = Check_software(logger, "vcftools1")
    beagle   = Check_software(logger, "beagle.jar")
    plink    = Check_software(logger, "plink")

    filterSNP   = FilterSNP(gatk, vcftools, beagle, plink, args, logger)
    #filterSNP.gatkfilter()
    filterSNP.missingfilter()
    filterSNP.imputation()
    #filterSNP.result()
    filterSNP.bialleleMafilter()
    filterSNP.linkagedisequilibrium()

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:
python3 {scriptpath}/FilterSNP.py <args> <args>....
NOTE:
This script was used to filterout SNPs that failed to pass GATK criterion and with missing rate in population above threshold.
'''.format(scriptpath = scriptpath))

    parse.add_argument('-invcf', '--inputvcf', required = True, dest = "invcf", help = "raw genotyping vcf", metavar = "Raw genotyping vcf file", type = str, nargs = '?')
    parse.add_argument('-logfile', '--logfile', required = True, dest = "logfile", help = "Log file to record procedures of processing of this script", metavar = "Log file to record procedures of processing of this script", type = str, nargs = '?')
    parse.add_argument('-ref', '--refgenome', required = True, dest = "refgenome", help = "ref genome in fa format", type = str, nargs = '?')
    parse.add_argument('-rate', '--threshold', required = False, dest = "missingthreshold", help = "SNPs max missing threshold in population", metavar = "SNPs max missing threshold in population", type = float, default = 0.01, nargs = '?')
    parse.add_argument('-tempdir', '--tempdir', required = True, dest = "tempdir", help = "Temp dir for all java programs", type = str, nargs = '?')
    parse.add_argument('-nt', '--nthreads', required = False, dest = "nthreads", help = "Number of threads for one run", type = int, default = 4, nargs = '?')
    parse.add_argument('-redir', '--resultdir', required = True, dest = "resultdir", help = "Result for stat files", type = str, nargs = '?')
    parse.add_argument('-maf', '--maf', required = False, dest = "maf", help = "Minimun allele frequency", type = float, default = 0.01, nargs = '?')
    parse.add_argument('-ld', '--LD', required = False, dest = "LD", help = "Threshold to define LD", type = float, default = 0.0001, nargs = '?')

    args = parse.parse_args()

    main(args)

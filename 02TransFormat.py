#!/usr/bin/env python3

import re,time,os,sys
import traceback
import argparse
import json
import subprocess
import configparser
from functools import reduce
import logging

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

class TransFormat():
    def __init__(self, args, logger):
        self.snpfile = args.insnpvcf
        self.outsnp  = args.outsnp
        self.logger  = logger
        self.genotypedict = {"0/0":"A","0/1":"H","1/1":"B","./.":"-"}

        self.snpdir = os.path.dirname(os.path.abspath(self.snpfile))

    def transformat(self):
        self.logger.info("===== Transformatting is starting =====")
        transfile = self.snpdir + "/snp.imputation.trans.txt"
        o = open(transfile, 'w+')
        with open(self.snpfile, 'r') as snpfile:
            for line in snpfile:
                translist = list()
                if not line.startswith("##"):
                    contentlist = line.strip().split("\t")
                    if line.startswith("#CHROM"):
                        o.write("MarkerID\tCHROM" + "\t" + contentlist[1] + "\t" + "\t".join(contentlist[9:]) + "\n")
                    else:
                        chrom = contentlist[0]
                        pos   = contentlist[1]
                        samplenum = len(contentlist[9:])
                        for i in range(samplenum):
                            genotype = contentlist[9 + i].split(":")[0]
                            translist.append(self.genotypedict[genotype])
                        MarkerID = chrom + ":" + pos
                        o.write(MarkerID + "\t" + chrom + "\t" + pos + "\t" + "\t".join(translist) + "\n")
        o.close()

    def filterhete(self):
        transfile = self.snpdir + "/snp.imputation.trans.txt"
        o = open(self.outsnp, 'w+')
        with open(transfile, 'r') as trans:
            for line in trans:
                if line.startswith("MarkerID"):
                    o.write(line)
                else:
                    contentlist = line.strip().split("\t")
                    patenal = contentlist[3]
                    matenal = contentlist[4]
                    if patenal != matenal:
                        if patenal != "H" and matenal != "H":
                            o.write(line)
        o.close()

def main(args):
    vcffile  = args.insnpvcf
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

    transformat   = TransFormat(args, logger)
    transformat.transformat()
    transformat.filterhete()

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:
python3 {scriptpath}/02TransFormat.py <args> <args>....

NOTE:
This script was used to convert genotypes to specific chracters in order to downstream analysis.

Normally, transformation is as follows:

0/0 =>  A,
1/1 =>  B,
0/1 =>  H.

Parental genotypes should be homozygous theoretically, so the SNP of materal or paternal with H will be filtered.
'''.format(scriptpath = scriptpath))

    parse.add_argument('-insnpvcf', '--inputvcf', required = True, dest = "insnpvcf", help = "snp vcf imputated and noise filtered", metavar = "vcf imputated and noise filtered", type = str, nargs = '?')
    parse.add_argument('-outsnp', '--outputsnp', required = True, dest = "outsnp", help = "snp genotype file with A,B,H converted", metavar = "snp genotype file with A,B,H converted", type = str, nargs = '?')
    parse.add_argument('-logfile', '--logfile', required = True, dest = "logfile", help = "Log file to record procedures of processing of this script", metavar = "Log file to record procedures of processing of this script", type = str, nargs = '?')

    args = parse.parse_args()

    main(args)

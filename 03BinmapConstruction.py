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

class BinConstruct():
    def __init__(self, args, logger):
        self.transfile = args.transfile
        self.binsize   = args.binsize
        self.binnorm   = args.binnorm
        self.logger    = logger
        self.outfile   = args.outfile

    def binconstruct(self):
        o = open(self.outfile, 'w+')
        test = open("test.bin.txt", "w+")
        o.write("BINMARKER\tCHR\tSTART\tEND\tSNPS\tGENOTYPES\n")
        infodict = dict()
        chromdict = dict()
        rowindex = 0
        title = ""
        with open(self.transfile, 'r') as f:
            title = f.readline()
            for line in f:
                infodict[rowindex] = dict()
                contentlist = line.strip().split("\t")
                chrom   = contentlist[1]
                if chrom not in chromdict:
                    chromrowindex = 0
                    chromdict[chrom] = dict()
                    chromdict[chrom][chromrowindex] = line
                else:
                    chromrowindex += 1
                    chromdict[chrom][chromrowindex] = line
                colnums = len(contentlist)
                for i in range(colnums):
                    infodict[rowindex][i] = contentlist[i]
                rowindex += 1
        for chrom in chromdict:
            snpnums = len(chromdict[chrom])
            slindict = dict()
            if snpnums < self.binsize:
                pass
            else:
                residue = (snpnums - self.binsize)%2
                #titlelist = title.strip().split("\t")
                #colnums = len(titlelist)
                #for i in range(5, colnums):

                if residue == 0:
                    for i in range(0, snpnums - self.binsize, 2):
                        sliname   = chrom + "-" + "Marker" + str(i)
                        slinstart = chromdict[chrom][i].strip().split("\t")[2]
                        slinend   = chromdict[chrom][i+self.binsize].strip().split("\t")[2]
                        colnums   = len(chromdict[chrom][i].strip().split("\t"))
                        slinsnps  = list()
                        slingenotypes = list()
                        genotypes = list()
                        for j in range(5, colnums):
                            for k in range(self.binsize):
                                contents = chromdict[chrom][i+k].strip().split("\t")
                                snpname  = contents[0]
                                genotype = contents[j]
                                slinsnps.append(snpname)
                                genotypes.append(genotype)
                            Anums = genotypes.count("A")
                            Bnums = genotypes.count("B")
                            Hnums = genotypes.count("H")
                            if Anums >= self.binnorm:
                                slingenotypes.append("A")
                            elif Bnums >= self.binnorm:
                                slingenotypes.append("B")
                            else:
                                slingenotypes.append("H")
                        slindict[sliname] = list()
                        slindict[sliname].extend([sliname, chrom, slinstart, slinend, slinsnps, slingenotypes])
                else:
                    for i in range(0, snpnums - self.binsize - 1, 2):
                        sliname   = chrom + "-" + "Marker" + str(i)
                        slinstart = chromdict[chrom][i].strip().split("\t")[2]
                        colnums   = len(chromdict[chrom][i].strip().split("\t"))
                        slinsnps  = list()
                        slingenotypes = list()
                        genotypes = list()
                        if i == snpnums - self.binsize - 2:
                            slinend = chromdict[chrom][i+self.binsize+1].strip().split("\t")[2]
                            for j in range(5, colnums):
                                for k in range(self.binsize + 1):
                                    contents = chromdict[chrom][i+k].strip().split("\t")
                                    snpname  = contents[0]
                                    genotype = contents[j]
                                    slinsnps.append(snpname)
                                    genotypes.append(genotype)
                                Anums = genotypes.count("A")
                                Bnums = genotypes.count("B")
                                Hnums = genotypes.count("H")
                                if Anums >= self.binnorm:
                                    slingenotypes.append("A")
                                elif Bnums >= self.binnorm:
                                    slingenotypes.append("B")
                                else:
                                    slingenotypes.append("H")
                            slindict[sliname] = list()
                            slindict[sliname].extend([sliname, chrom, slinstart, slinend, slinsnps, slingenotypes])
                        else:
                            slinend = chromdict[chrom][i+self.binsize].strip().split("\t")[2]
                            for j in range(5, colnums):
                                for k in range(self.binsize):
                                    contents = chromdict[chrom][i+k].strip().split("\t")
                                    snpname  = contents[0]
                                    genotype = contents[j]
                                    slinsnps.append(snpname)
                                    genotypes.append(genotype)
                                Anums = genotypes.count("A")
                                Bnums = genotypes.count("B")
                                Hnums = genotypes.count("H")
                                if Anums >= self.binnorm:
                                    slingenotypes.append("A")
                                elif Bnums >= self.binnorm:
                                    slingenotypes.append("B")
                                else:
                                    slingenotypes.append("H")
                            slindict[sliname] = list()
                            slindict[sliname].extend([sliname, chrom, slinstart, slinend, slinsnps, slingenotypes])
            slinnums = len(slindict)
            for sliname in slindict:
                test.write(slindict[sliname][0] + "\t" + slindict[sliname][1] + "\t" + slindict[sliname][2] + "\t" + slindict[sliname][3] + "\t" + "\t".join(slindict[sliname][4]) + "\t" + "\t".join(slindict[sliname][5]) + "\n" )
            tmpkey = chrom + "-" + "Marker0"
            binnum = 0
            bindict = dict()
            for i in range(1, slinnums):
                key1 = chrom + "-" + "Marker" + str(i*2)
                if slindict[key1][5] == slindict[tmpkey][5]:
                    #slindict[tmpkey] = list()
                    sliname = slindict[key1][0]
                    chrom   = slindict[key1][1]
                    slinstart = slindict[tmpkey][2]
                    slinend = slindict[key1][3]
                    slinsnps = slindict[tmpkey][4] + slindict[key1][4]
                    slingenotypes = slindict[key1][5]
                    tmpkey = key1
                    slindict[tmpkey] = list()
                    slindict[tmpkey].extend([sliname, chrom, slinstart, slinend, slinsnps, slingenotypes])
                else:
                    binmark = chrom + "-" + "Bin" + str(binnum)
                    chrom = slindict[tmpkey][1]
                    binstart = slindict[tmpkey][2]
                    binend   = slindict[key1][2]
                    binsnps  = slindict[tmpkey][4]
                    bingenotypes = slindict[tmpkey][5]
                    bindict[binmark] = list()
                    bindict[binmark].extend([binmark, chrom, binstart, binend, binsnps, bingenotypes])
                    binnum += 1
                    linecontent = bindict[binmark][0] + "\t" + bindict[binmark][1] + "\t" + str(bindict[binmark][2]) + "\t" + str(bindict[binmark][3]) + "\t" + ",".join(bindict[binmark][4]) + "\t" + ",".join(bindict[binmark][5]) + "\n"
                    o.write(linecontent)
            #for binmark in bindict:
                #binlength = int(bindict[binmark][3]) - int(bindict[binmark][2])
                #if binlength >= 300:
                    #linecontent = bindict[binmark][0] + "\t" + bindict[binmark][1] + "\t" + str(bindict[binmark][2]) + "\t" + str(bindict[binmark][3]) + "\t" + ",".join(bindict[binmark][4]) + "\t" + ",".join(bindict[binmark][5]) + "\n"
                    #o.write(linecontent)
        o.close()

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

    binconstruct = BinConstruct(args, logger)
    binconstruct.binconstruct()

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:
python3 {scriptpath}/03BinmapConstruction.py <args> <args>....

NOTE:
This script was used to reduce snps adjacent to each other with same genotypes into bins.

'''.format(scriptpath = scriptpath))

    parse.add_argument('-transfile', '--transfile', required = True, dest = "transfile", help = "file with genotypes in converted characters", metavar = "file with genotypes in converted characters", type = str, nargs = '?')
    parse.add_argument('-outfile', '--outfile', required = True, dest = "outfile", help = "outfile with bin marker info", metavar = "outfile with bin marker info", type = str, nargs = '?')
    parse.add_argument('-binsize', '--binsize', required = True, dest = "binsize", help = "bin size to define a slinding window", metavar = "bin size to define a slinding window", type = int, nargs = '?')
    parse.add_argument('-binnorm', '--binorm', required = True, dest = "binnorm", help = "strandard or threshold to define genotype in a window for one individual", type = int, nargs = '?')
    parse.add_argument('-logfile', '--logfile', required = True, dest = "logfile", help = "Log file to record procedures of processing of this script", metavar = "Log file to record procedures of processing of this script", type = str, nargs = '?')

    args = parse.parse_args()

    main(args)

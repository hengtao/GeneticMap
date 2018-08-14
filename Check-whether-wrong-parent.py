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

def CheckWrongParent(snpfile, outfile1, outfile2):
	genotypedict = {"0/0":"A","0/1":"H","1/1":"B"}
	o1 = open(outfile1, 'r')
	o2 = open(outfile2, 'r')
	
	with open(snpfile, 'r') as f:
		for line in f:
			if line.startswith("#"):
				pass
			else:
				(NN, AA, BB, HH, EE) = (0,0,0,0,0)
				contentlist = line.strip().split("\t")
				alt = contentlist[4]
				altnum = len(alt)
				if altnum > 1:
					par = contentlist[9].split(":")[0]
					mar = contentlist[10].split(":")[0]
					progencynum = len(contentlist) - 11
					if mar != "./." and mar not in genotypedict and par != "./." and par in genotypedict:
						for i in range(progencynum):
							geno = contentlist[11+i].split(":")[0]
							if geno == "./.":
								NN += 1
							elif geno == "0/0":
								AA += 1
							elif geno == "0/1":
								HH += 1
							elif geno == "1/1":
								BB += 1
							else:
								EE += 1
						if NN != progencynum:
							numlist = [AA,HH,BB,EE,NN]
							o1.write("\t".join([str(num) for num in numlist]) + "\t" + "\t".join(contentlist) + "\n")
					elif mar != "./." and mar in genotypedict and par != "./." and par not in genotypedict:
						for i in range(progencynum):
							geno = contentlist[11+i].split(":")[0]
							if geno == "./.":
								NN += 1
							elif geno == "0/0":
								AA += 1
							elif geno == "0/1":
								HH += 1
							elif geno == "1/1":
								BB += 1
							else:
								EE += 1
						if NN != progencynum:
							numlist = [AA,HH,BB,EE,NN]
							o2.write("\t".join([str(num) for num in numlist]) + "\t" + "\t".join(contentlist) + "\n")
	o1.close()
	o2.close()
							
				
def main(args):
    logfile  = args.logfile
	snpfile  = args.invcf
	outfile1 = args.parout
	outfile2 = args.marout

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
	
	CheckWrongParent(snpfile, outfile1, outfile2)
	
if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:
python3 {scriptpath}/Check-whether-wrong-parent.py <args> <args>....
NOTE:
This script was used to check which parent was wrongly chose for GBS program.
'''.format(scriptpath = scriptpath))
 
	parse.add_argument('-snpvcf', '--insnpvcf', required = True, dest = "invcf", help = "Filtered genotyping vcf", metavar = "Only contain snps that pass gatk filteration standard", type = str, nargs = '?')
    parse.add_argument('-parout', '--parout', required = True, dest = "parout", help = "Evidence supporting the correct male parent", metavar = "Multi-alt snp, geno of male parent is 0/0, 0/1 or 1/1, while female is 0/2, 1/2 or 2/2", type = str, nargs = '?')
	parse.add_argument('-marout', '--marout', required = True, dest = "marout", help = "Evidence supporting the correct female parent", metavar = "Multi-alt snp, geno of female parent is 0/0, 0/1 or 1/1, while male is 0/2, 1/2 or 2/2", type = str, nargs = '?')

	args = parse.parse_args()

    main(args)
	
	

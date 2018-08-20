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
import http.cookiejar
import urllib.request
import urllib.parse
import urllib2
from bs4 import BeautifulSoap as bs

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

def EvidenceLevel(evidence, queryvcf, outputfile):
	FuncrefGene = re.compile(r"Func.refGene=(.*?);")
	GenerefGene = re.compile(r"Gene.refGene=(.*?);")
	ExonicFuncrefGene = re.compile(r"ExonicFunc.refGene=(.*?);")
	AAChangerefGene   = re.compile(r"AAChange.refGene=(.*?);")
	condidategenes = list()
	o = open(outputfile, 'w+')
	with open(queryvcf) as f:
		for line in f:
			contentlist = line.strip().split("\t")	
			chrom = contentlist[0]
			posit = contentlist[1]
			refer = contentlist[3]
			alter = contentlist[4]
			infor = contentlist[7]
			MatchFuncrefGene = re.findall(FuncrefGene, infor)[1]
			MatchGenerefGene = re.findall(GenerefGene, infor)[1]
			MatchExonicFuncrefGene = re.findall(ExonicFuncrefGene, infor)[1]
			MatchAAChangerefGene   = re.findall(AAChangerefGene, infor)[1]
			if MatchFuncrefGene in evidence:
				o.write(MatchGenerefGene + "\t" + chrom + "\t" + posit + "\t" + refer + "\t" + alter + "\t" + MatchFuncrefGene + "\t" + MatchExonicFuncrefGene + "\t" +  MatchAAChangerefGene + "\n")
			if MatchFuncrefGene == "Exonic":
				if MatchExonicFuncrefGene == "nonsynonymous_SNV" or MatchExonicFuncrefGene.startswith("stop") or MatchExonicFuncrefGene.startswith("frameshift"):
					if MatchGenerefGene not in condidategenes:
						condidategenes.append(MatchGenerefGene)
	o.close()
	return MatchGenerefGene
						
def VersionConvert(args):
	versionFile = args.versionFile
	evidence = args.evidence
	queryvcf = args.queryvcf
	outputfile = args.outputfile
	versoionname = args.versionName
	genenamefile = args.genenamefile
	querygene  = EvidenceLevel(evidence, queryvcf, outputfile)
	condidategenes = list()
	targetgenes = list()
	o = open(querygene, 'w+')
	with open(versionFile, 'r') as f:
		for line in f:
			for gene in querygene:
				if re.findall(gene, line):
					del(querygene[0])
					linecontents = line.strip().split("\t")[1:]
					lineversionname  = [content.split("|")[1] if content.startswith(versoionname) for content in linecontents]
					linewithoutname  = [content.split("|")[1] if not content.startswith(versoionname) for content in linecontents]	
					if gene in lineversionname:
						condidategenes.extend(linewithoutname)
						o.write(gene + "\t" + ",".join(linewithoutname))
	o.close()
	for gene in condidategenes:
		if gene not in targetgenes:
			targetgenes.append(gene)
	return targetgenes
	
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

	evidence = args.evidence
	queryvcf = args.queryvcf
	outputfile = args.outputfile
	
    if args.versionFile:
		targetgenes = VersionConvert(args)
	else:
		targetgenes = EvidenceLevel(evidence, queryvcf, outputfile)

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:
python3 {scriptpath}/FilterSNP.py <args> <args>....
NOTE:
This script was used to filterout SNPs that failed to pass GATK criterion and with missing rate in population above threshold.
'''.format(scriptpath = scriptpath))

    parse.add_argument('-query', '--queryvcf', required = True, dest = "queryvcf", help = "vcf file without intergenic regions info", metavar = "vcf file without intergenic regions info", type = str, nargs = '?')
    parse.add_argument('-outfile', '--out', required = True, dest = "outputfile", help = "outputfile record condidates SNP/Indel", metavar = "outputfile record condidates SNP/Indel", type = str, nargs = '?')
    parse.add_argument('-genenamefile', '--namefile', required = True, dest = "genefile", help = "condidate genes with their corresponding names of another version", type = str, nargs = '?')

	parse.add_argument('-logfile', '--logfile', required = True, dest = "logfile", help = "Log file to record procedures of processing of this script", metavar = "Log file to record procedures of processing of this script", type = str, nargs = '?')
	parse.add_argument('-evidence', '--evidence', required = True, dest = "evidence", default = "exonic downstream intronic splicing upstream UTR3 UTR5", 
					   help = "Positon where SNP/indel located in will be further studied", type = str, nargs = '*')
	versiongroup = parse.add_argument_group("verion conversion")
	versiongroup.add_argument('-versionFile', '--versionFile', required = False, dest = "versionFile", help = "Responding relationship between two versions", type = str, nargs = '?')
	versiongroup.add_argument('-versionname', '--versionname', required = False, dest = "versionName", help = "version name should be consistent with that in verionfile", type = str, nargs = '?')
	
	args = parse.parse_args()

    main(args)

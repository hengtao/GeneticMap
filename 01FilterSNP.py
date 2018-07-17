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

class FilterSNP():
    def __init__(self, gatk, vcftools, beagle, plink,  args, logger):
        #TOOLS
        self.gatk = gatk
        self.vcftools = vcftools
        self.beagle = beagle
        self.plink  = plink
        #INPUT&&OUTPUT
        self.invcf = args.invcf
        #self.outvcf = args.outvcf
        self.logger = logger
        self.threshold = args.missingthreshold
        self.tempdir = args.tempdir
        self.refgenome = args.refgenome
        self.nthreads = args.nthreads
        self.maf = args.maf
        self.LD = args.LD

        self.resultdir = args.resultdir
        self.vcfdir = os.path.dirname(os.path.abspath(self.invcf))
        self.scriptpath = os.path.split(os.path.realpath(__file__))[0]

    def gatkfilter(self):
        samplesline = os.popen("""grep "#CHROM" {invcf}""".format(invcf = self.invcf)).read()
        samplesnum  = len(samplesline.strip().split("\t")) - 9
        print(samplesline)
        outsnpvcf   = self.vcfdir + "/snp.raw.vcf"
        outindelvcf = self.vcfdir + "/indel.raw.vcf"
        outsnpfilter = self.vcfdir + "/snp.gatkfilter.vcf"
        outindelfilter = self.vcfdir + "/indel.gatkfilter.vcf"
        outsnpnomissing = self.vcfdir + "/snp.missfilter.vcf"

        cmd = "java -XX:ParallelGCThreads=4 -Djava.io.tmpdir={tempdir} -Xmx10g -jar {gatk} -T SelectVariants -R {refgenome} -V {invcf} -selectType SNP -o {outsnpvcf} \n".format(tempdir = self.tempdir, gatk = self.gatk, refgenome = self.refgenome, invcf = self.invcf, outsnpvcf = outsnpvcf)
        cmd += "java -XX:ParallelGCThreads=4 -Djava.io.tmpdir={tempdir} -Xmx10g -jar {gatk} -T SelectVariants -R {refgenome} -V {invcf} -selectType INDEL -o {outindelvcf} \n".format(tempdir = self.tempdir, gatk = self.gatk, refgenome = self.refgenome, invcf = self.invcf, outindelvcf = outindelvcf)
        cmd += "java -XX:ParallelGCThreads=4 -Djava.io.tmpdir={tempdir} -Xmx30g -jar {gatk} -T VariantFiltration -R {refgenome} -V {outsnpvcf} --filterExpression \" QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 4.0 || ReadPosRankSum < -8.0 \" --missingValuesInExpressionsShouldEvaluateAsFailing  --clusterWindowSize 10 --clusterSize 4  --filterName  snp_filter  -o {outsnpfilter} \n ".format(tempdir = self.tempdir, gatk = self.gatk, refgenome = self.refgenome, outsnpfilter = outsnpfilter, outsnpvcf = outsnpvcf)
        cmd += "java -XX:ParallelGCThreads=4 -Djava.io.tmpdir={tempdir} -Xmx30g -jar {gatk} -T VariantFiltration -R {refgenome} -V {outindelvcf} ".format(tempdir = self.tempdir, gatk = self.gatk, refgenome = self.refgenome, outindelvcf = outindelvcf)
        if int(samplesnum)> 10:
            cmd += "--filterExpression \" QUAL < 30.0 || QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || MQ < 40.0 || MQRankSum < -12.5 || InbreedingCoeff < -0.8\" "
        else:
            cmd += "--filterExpression \" QUAL < 30.0 || QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || MQ < 40.0 || MQRankSum < -12.5 \" "
        cmd += "--missingValuesInExpressionsShouldEvaluateAsFailing  --filterName  indel_filter  -o {outindelfilter} ".format(outindelfilter = outindelfilter)
        Cmd_and_Time(cmd, self.logger)

    def missingfilter(self):
        outsnpfilter = self.vcfdir + "/snp.gatkfilter.vcf"
        outsnpnomissing = self.vcfdir + "/snp.missfilter"
        cmd = "{vcftools} --vcf {outsnpfilter} --max-missing {threshold} --recode --recode-INFO-all --out {outsnpnomissing}".format(vcftools = self.vcftools, outsnpfilter = outsnpfilter, threshold = self.threshold, outsnpnomissing = outsnpnomissing)
        Cmd_and_Time(cmd, self.logger)

    def imputation(self):
        outsnpnomissing = self.vcfdir + "/snp.missfilter.recode.vcf"
        snpimpute = self.vcfdir + "/snp.impute"
        indelimpute = self.vcfdir + "/indel.impute"
        outindelfilter = self.vcfdir + "/indel.gatkfilter.vcf"
        headerfile = self.vcfdir + "/vcfheader.txt"
        snpnoheader = self.vcfdir + "/snp.impute.noheader.vcf"
        indnoheader = self.vcfdir + "/indel.impute.noheader.vcf"
        mergeimpute = self.resultdir + "/merge.impute.vcf"
        mergesort   = self.resultdir + "/merge.impute.sort.vcf"
        cmd = "java -XX:ParallelGCThreads=4 -Djava.io.tmpdir={tempdir} -Xmx30g -jar  {beagle} gtgl={outsnpnomissing} out={snpimpute} nthreads={nthreads} impute=true\n".format(tempdir = self.tempdir, beagle = self.beagle, outsnpnomissing = outsnpnomissing, snpimpute = snpimpute, nthreads = self.nthreads)
        cmd += "java -XX:ParallelGCThreads=4 -Djava.io.tmpdir={tempdir} -Xmx20g -jar {beagle} gtgl={outindelfilter} out={indelimpute} nthreads={nthreads} impute=true".format(tempdir = self.tempdir, beagle = self.beagle, outindelfilter = outindelfilter, indelimpute = indelimpute, nthreads = self.nthreads)
        Cmd_and_Time(cmd, self.logger)
        if not os.path.exists(headerfile):
            os.system("zcat {snpimpute}.vcf.gz |grep  '^#' > {headerfile}".format(snpimpute = snpimpute, headerfile = headerfile))
        os.system("zcat {snpimpute}.vcf.gz |grep -v '^#' > {snpnoheader}".format(snpimpute = snpimpute, snpnoheader = snpnoheader))
        os.system("zcat {indelimpute}.vcf.gz |grep -v '^#' > {indnoheader}".format(indelimpute = indelimpute, indnoheader = indnoheader))
        os.system("cat {snpnoheader} {indnoheader} | sort -nk 1 > {redir}/temp.txt".format(snpnoheader = snpnoheader, indnoheader = indnoheader, redir = self.resultdir))
        os.system("cat {header} {redir}/temp.txt > {mergeimpute}".format(header = headerfile, redir = self.resultdir, mergeimpute = mergeimpute))
        os.system("python3 {Bin}/../SNPScripts/sort_imputation_vcf.py -i {mergeimpute} -o {mergeimputesort}".format(Bin = self.scriptpath, mergeimpute = mergeimpute, mergeimputesort = mergesort))
        os.system("rm -f {header} {snpnoheader} {indelnoheader} {redir}/temp.txt {mergeimpute} {snpimpute}.log {indelimpute}.log".format(header = headerfile, snpnoheader = snpnoheader, indelnoheader = indnoheader, redir = self.resultdir, mergeimpute = mergeimpute, snpimpute = snpimpute, indelimpute = indelimpute))

    def result(self):
        self.logger.info("===== Extract VCF information is starting =====")
        outsnpfilter = self.vcfdir + "/snp.gatkfilter.vcf"
        outindelfilter = self.vcfdir + "/indel.gatkfilter.vcf"

        cmd = "cd {RESdir} && {Bin}/../bin/ExtractVCF.py SNP {outsnpfilter} \n".format(RESdir = self.resultdir, Bin = self.scriptpath, outsnpfilter = outsnpfilter)
        cmd += "cd {RESdir} && {Bin}/../bin/ExtractVCF.py INDEL {outindelfilter} \n".format(RESdir = self.resultdir, Bin = self.scriptpath, outindelfilter = outindelfilter)
        cmd += "cd {RESdir} && {vcftools} --vcf snp.filter.vcf --SNPdensity 10000 --out 10k \n".format(RESdir = self.resultdir, vcftools = self.vcftools)
        cmd += "cd {RESdir} && perl {Bin}/../bin/IndelStat.pl indel.list 20 > indel.length.stat ".format(RESdir = self.resultdir, Bin = self.scriptpath)
        Cmd_and_Time(cmd, self.logger)

    def bialleleMafilter(self):
        self.logger.info("===== Filterout SNPs without Bi-alleles is starting =====")
        snpimpute = self.vcfdir + "/snp.impute.vcf.gz"
        snpbiallele = self.vcfdir + "/snp.impute.biallele"
        cmd = "{vcftools} --gzvcf {snpimpute} --min-alleles 2 --max-alleles 2 --maf {maf}  --recode --recode-INFO-all --out {snpbiallele}".format(vcftools = self.vcftools, snpimpute = snpimpute, maf = self.maf, snpbiallele = snpbiallele)
        Cmd_and_Time(cmd, self.logger)

    def linkagedisequilibrium(self):
        self.logger.info("===== Filterout SNPs with LD is starting =====")
        snpbiallele = self.vcfdir + "/snp.impute.biallele.recode.vcf"
        snpldvcf  = self.vcfdir + "/snp.impute.LD.vcf"
        cmd = "cd {vcfdir} && {vcftools} --vcf {snpbiallele} --plink --out LD \n".format(vcfdir = self.vcfdir, vcftools = self.vcftools, snpbiallele = snpbiallele)
        cmd += "cd {vcfdir} && {plink} --file LD --make-bed --out gbs --noweb --geno 0.05 --maf {maf} --hwe {LD}".format(vcfdir = self.vcfdir, plink = self.plink, maf = self.maf, LD = self.LD)
        Cmd_and_Time(cmd, self.logger)
        chrposlist = list()
        with open(self.vcfdir + "/gbs.bim", 'r') as bim:
            for line in bim:
                chrpos = line.strip().split("\t")[1]
                if chrpos not in chrposlist:
                    chrposlist.append(chrpos)
        o = open(snpldvcf, 'w+')
        with open(snpbiallele ,'r') as f:
            for line in f:
                if line.startswith("#"):
                    o.write(line)
                else:
                    contentlist = line.strip().split("\t")
                    chrom = contentlist[0]
                    posit = contentlist[1]
                    chrpos = chrom + ":" + posit
                    if chrpos in chrposlist:
                        o.write(line)

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


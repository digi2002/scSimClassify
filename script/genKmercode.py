import sys
import multiprocessing
from multiprocessing import Pool
import numpy as np
import glob
import os
import subprocess

# generate scaled features

def checkhashfile(nbit, zeroper, topper, codefilelist, abunfilelist):
    if len(codefilelist) != len(abunfilelist):
        print('Not all the abunfiles have the corresponding codefiles')
        return False
    for hashfile in codefilelist:
        hashfile = hashfile.replace('.hashcode' + str(nbit) + '_' + str(zeroper) + '_' + str(topper), '.abun.npz')
        if hashfile not in set(abunfilelist):
            print(hashfile + ' does not have the conresponding abunfile')
            return False
    return True


def selCode(nbit, zeroper, topper, hashdir, abunfilelist, outcodefile, outcodefile10):
    codefilelist = []
    codeCountdict = {}
    selCodeset = set()
    codeset = set()
    for filename in glob.iglob(hashdir + '*.hashcode' + str(nbit) + '_' + str(zeroper) + '_' + str(topper),
                               recursive=True):  # attention
        hashfile = filename.replace(hashdir, '')
        codefilelist.append(hashfile)
    if checkhashfile(nbit, zeroper, topper, codefilelist, abunfilelist) == False:
        print('the number of codefiles is not as the same as abundence files')
        return False
    for hashfile in codefilelist:
        f = open(hashdir + hashfile)
        for line in f:
            kmerid, code = line.strip().split('\t')
            if code not in codeCountdict.keys():
                codeCountdict[code] = 0
            codeCountdict[code] = codeCountdict[code] + 1
    for code, cnt in codeCountdict.items():
        codeset.add(code)
        if cnt>5:
            selCodeset.add(code)
    # output sorted codefile
    sortedSelcode = sorted(selCodeset, key=int)
    sortedCode = sorted(codeset, key=int)
    outcode10 = open(outcodefile10, 'w')
    outcode = open(outcodefile, 'w')
    for code in sortedSelcode:
        outcode10.write(str(code) + '\n')
    for code in sortedCode:
        outcode.write(str(code) + '\n')
    return set(sortedSelcode)

def genKmercodedict(nbit, zeroper, topper, hashdir, selCodeset, outfile10, outfile):
    print('selcodeset: ' + str(len(selCodeset)))
    output10 = open(outfile10, 'w')
    output = open(outfile, 'w')
    codefilelist = []
    for filename in glob.iglob(hashdir + '*.hashcode' + str(nbit) + '_' + str(zeroper) + '_' + str(topper),
                               recursive=True):  # attention
        codefilelist.append(filename)
    for filename in codefilelist:
        f = open(filename)
        for line in f:
            kmerid, code = line.strip().split('\t')
            output.write(line.strip() + '\n')
            if code in selCodeset:
                output10.write(line.strip() + '\n')







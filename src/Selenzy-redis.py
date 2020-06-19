#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 10:53:03 2017

@author: jerrywzy, Pablo Carbonell
"""
import re
from os import path as os_path
from os import system as os_system
from os import rename as os_rename
from json import load as json_load
from csv import reader as csv_reader
from csv import writer as csv_writer
from csv import QUOTE_ALL as csv_QUOTE_ALL
from numpy import zeros as np_zeros
from pandas import read_csv as pd_read_csv
from pandas import Series as pd_Series
from rdkit.Chem import AllChem, Draw
import quickRsim
from CRedisDict import CRedisDict, wait_for_redis
from redis import StrictRedis
from brs_utils import print_OK, print_FAILED
from time import time

class preLoad(object):
    """ Container of precomputed data """
    def __init__(self, datadir, sm='file'):
        self.datadir = datadir

        self.attr_names = {
            'fasta': ['sequence', 'names', 'descriptions', 'fulldescriptions', 'osource', 'pexistence'],
            'fp': ['fp']
            }
        for attr_name in self.attr_names:
            setattr(self, attr_name, {})

        self.store_mode = sm
        if self.store_mode!='file':
            self.redis = StrictRedis(host=self.store_mode, port=6379, db=0, decode_responses=True)
            if not wait_for_redis(self.redis, 30):
                self.logger.critical("Database "+self.store_mode+" is not reachable")
                exit()
            # Re-init attributes
            for cat in self.attr_names:
                for attr_name in self.attr_names[cat]:
                    setattr(self, attr_name, CRedisDict(attr_name, self.redis))

        self.MnxToUprot = {}
        self.UprotToMnx = {}

        print ("Acquiring databases...")
        self.readData()

    def readData(self):
        """ Read all data into memory """
        print('Get fasta...', end = '', flush=True)
        start = time()
        self._fasta('seqs.fasta')
        print_OK(time()-start)

        print('Get fpData...', end = '', flush=True)
        start = time()
        self._fpData()
        print_OK(time()-start)

        print('Get seqData...', end = '', flush=True)
        start = time()
        self._seqData(['reac_seqs.tsv', 'upclst.json', 'clstrep.json', "seq_org.tsv", "org_lineage.csv"])
        print_OK(time()-start)

        print('Get reacData...', end = '', flush=True)
        start = time()
        self._reacData(['reac_smi.csv','reac_xref.tsv',"rxn_consensus_20160612.txt",
                        'reac_prop.tsv', 'brenda-mnxref2.tsv', 'sabiork-mnxref2.tsv', 'reac_xref_v3.tsv'])
        print_OK(time()-start)

    def _fasta(self, ffile="seqs.fasta"):
        if self.store_mode=='file':
            (self.sequence,
             self.names,
             self.descriptions,
             self.fulldescriptions,
             self.osource,
             self.pexistence) = self.readFasta(ffile)
        else:
            all_set = True
            for attr_name in self.attr_names['fasta']:
                if not CRedisDict.exists(self.redis, attr_name):
                    all_set = False
                    break
            if not all_set:
                (sequence,
                 names,
                 descriptions,
                 fulldescriptions,
                 osource,
                 pexistence) = self.readFasta(ffile)
                 # Set attributes
                for attr_name in self.attr_names['fasta']:
                    if not CRedisDict.exists(self.redis, attr_name):
                        # If var is not a dict, then convert it
                        try:
                            setattr(self, attr_name, CRedisDict(attr_name, self.redis, locals()[attr_name]))
                        except TypeError:
                            setattr(self, attr_name, CRedisDict(attr_name, self.redis, {attr_name: locals()[attr_name]}))
            # self.sequence = CRedisDict('sequence', self.redis, sequence)
            # self.names = CRedisDict('names', self.redis, names)
            # self.descriptions = CRedisDict('descriptions', self.redis, descriptions)
            # self.fulldescriptions = CRedisDict('fulldescriptions', self.redis, fulldescriptions)
            # self.osource = CRedisDict('osource', self.redis, osource)
            # self.pexistence = CRedisDict('pexistence', self.redis, pexistence)

    def _fpData(self):
        do_it = True
        if self.store_mode!='file':
            if CRedisDict.exists(self.redis, 'fp'):
                do_it = False
        if do_it:
            for fpid in quickRsim.fingerprint():
                self.fp[fpid] = quickRsim.loadFingerprint(self.datadir, fpid)

    def _seqData(self, fl):
        with open(os_path.join(self.datadir, fl[0])) as f:
            for row in f:
                mnxr, db, seqid, source, ec = row.split('\t')
                if db == 'uniprot':
                    if mnxr not in self.MnxToUprot:
                        self.MnxToUprot[mnxr] = set()
                    self.MnxToUprot[mnxr].add(seqid)
                    if seqid not in self.UprotToMnx:
                        self.UprotToMnx[seqid] = set()
                    self.UprotToMnx[seqid].add(mnxr)

        with open(os_path.join(self.datadir, fl[1])) as f2:
            self.upclst = json_load(f2)

        with open(os_path.join(self.datadir, fl[2])) as f3:
            self.clstrep= json_load(f3)
        self.seqorg = self.seqOrganism(fl[3])
        self.tax = self.readTaxonomy(fl[4])

    def _reacData(self, smf):
        """ Transitional downgrade mapping for mnx v3 to v2.
        In the future everything will be based on v3.
        """
        smiFile = os_path.join(self.datadir, smf[0])
        rxnRefFile = os_path.join(self.datadir, smf[1])
        rxnConsensus = os_path.join(self.datadir, smf[2])
        rxnProp = os_path.join(self.datadir, smf[3])
        rxnBrenda = os_path.join(self.datadir, smf[4])
        rxnSabiork = os_path.join(self.datadir, smf[5])
        rxnv3 = None
        if len(smf) > 6:
            rxnv3 = os_path.join(self.datadir, smf[6])
        self.smir = {}
        if os_path.exists(smiFile):
            self.smir = reactionSmiles(smiFile)
        self.rxnref = {}
        if os_path.exists(rxnRefFile):
            self.rxnref = reactionXref(rxnRefFile, rxnBrenda, rxnSabiork,rxnv3)
        self.rxndir = {}
        if os_path.exists(rxnConsensus):
            self.rxndir = readRxnCons(rxnConsensus)
        self.ecrxn = {}
        self.rxnec = {}
        if os_path.exists(rxnProp):
            self.ecrxn, self.rxnec = readRxnProp(rxnProp)
        self.ecsmi = ecSmiles(self.ecrxn, self.smir, self.rxnref)

    def seqOrganism(self, fileSeqOrg):
        seqorg = {}
        for line in open(os_path.join(self.datadir, fileSeqOrg)):
            row = line.rstrip().split('\t')
            seqorg[ row[0] ] = (row[1], row[2])
        return seqorg

    def readTaxonomy(self, fileLineage):
        tax = {}
        with open(os_path.join(self.datadir, fileLineage)) as handler:
            for row in csv_reader(handler):
                tax[row[0]] = row
        return tax

    def readFasta(self, fileFasta, limit=None):

        from Bio import SeqIO

        sequence = {}
        descriptions = {}
        fulldescriptions = {}
        osource = {}
        pexistence = {}
        names = []
        seen = set()
        seen_add = seen.add

        for seq_record in SeqIO.parse(os_path.join(self.datadir, fileFasta), "fasta"):
            ming = seq_record.id
            try:
                idonly = re.search(r'\|(.*?)\|',ming)
                x = idonly.group(1)
            except:
                x = ming
            fulldesc = seq_record.description
            desc = fulldesc.rsplit('OS=')[0]
            try:
                shortdesc = " ".join(desc.split()[1:])
            except:
                shortdesc = desc
            try:
                orgsource = fulldesc.rsplit('OS=')[1]
            except:
                orgsource = '-'
            shortos = orgsource.rsplit('GN=')[0]
            if ',' in shortdesc:
                y = shortdesc.replace(",", ";")
            else:
                y = shortdesc
            pe = orgsource.rsplit('PE=')[1].split(' ')[0]
            if x not in seen:
                names.append(x)
                seen_add(x)
            myseq = seq_record.seq
            sequence[x]=(myseq)
            descriptions[x]= y
            osource[x]=shortos
            fulldescriptions[x] = fulldesc
            pexistence[x] = pe
            # A practical limit hardcoded
            if limit is not None:
                try:
                    if len(sequence) > limit:
                        break
                except:
                    continue

        return (sequence, names, descriptions, fulldescriptions, osource, pexistence)


def availableFingerprints():
    fp = quickRsim.fingerprint()
    return sorted(fp, key=lambda x: fp[x][4])

def sanitizeRxn(rxninfo, outrxn):
    """ It works both with the smiles string or a rxn file """
    try:
        if os_path.exists(rxninfo):
            if open(rxninfo).readline().startswith('$RXN'):
                rxn = AllChem.ReactionFromRxnFile(rxninfo)
            else:
                smarts =  open(rxninfo).readline().rstrip()
                try:
                    rxn = AllChem.ReactionFromSmarts(smarts, useSmiles=True)
                except:
                    rxn = AllChem.ReactionFromSmarts(smarts)
        smi = AllChem.ReactionToSmiles(rxn)
        with open(outrxn+'.smi', 'w') as handler:
            handler.write(smi)
        mdl = AllChem.ReactionToRxnBlock(rxn)
        with open(outrxn+'.rxn', 'w') as handler:
            handler.write(mdl)
        return smi
    except:
        return ''

def sanitizeSmarts(smarts, outrxn):
    """ It works both with the smiles string or a rxn file """
    try:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts, useSmiles=True)
        except:
            rxn = AllChem.ReactionFromSmarts(smarts)
        smi = AllChem.ReactionToSmiles(rxn)
        with open(outrxn+'.smi', 'w') as handler:
            handler.write(smi)
        mdl = AllChem.ReactionToRxnBlock(rxn)
        with open(outrxn+'.rxn', 'w') as handler:
            handler.write(mdl)
        return smi
    except:
        return ''

def display_reaction(rxninfo, outfolder, outname, marvin=False):
    """ It works both with the smiles string or a rxn file """
    if marvin:
        try:
            cmd = ['molconvert', 'svg:w500', rxninfo]
            job = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = job.communicate()
            outimage = os_path.join(outfolder, outname+'.svg')
            return outimage, (600, 400)
        except:
            return '', (0,0)
    else:
        try:
            if os_path.exists(rxninfo):
                if open(rxninfo).readline().startswith('$RXN'):
                    rxn = AllChem.ReactionFromRxnFile(rxninfo)
                else:
                    smarts =  open(rxninfo).readline()
                    rxn = AllChem.ReactionFromSmarts(smarts)
                outimage = os_path.join(outfolder, outname+'.png')
                im = Draw.ReactionToImage(rxn)
                size = im.size
                im.save(outimage)
                return outimage, size
            else:
                rxn = AllChem.ReactionFromSmarts(rxninfo)
                outimage = os_path.join(outfolder, outname+'.png')
                im = Draw.ReactionToImage(rxn)
                size = im.size
                im.save(outimage)
                return outimage, size
        except:
            return '', (0,0)

def taxDistance(tax, host, target):
    if host in tax and target in tax:
        hostLineage = set(tax[host])
        targetLineage = set(tax[target])
        distance = 1 + len(hostLineage ^ targetLineage)
        return distance
    else:
        return '-'

def seqScore(newscore=None):
    import string
    # Initial score
    vdict = {
        string.ascii_uppercase[9]: ('Reaction similarity:', 100.0, True),
        string.ascii_uppercase[8]: ('Sequence <a href="https://en.wikipedia.org/wiki/Conserved_sequence" target="_blank">conservation</a>:', 1.0, True),
        string.ascii_uppercase[4]: ('Sequence <a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank">taxonomic</a> distance:', -1.0, True),
        string.ascii_uppercase[7]: ('Uniprot <a href="http://www.uniprot.org/help/protein_existence" target="_blank">protein evidence</a>:', -0.1, True),
        string.ascii_uppercase[12]: ('Percentage helices:', 0.0, False),
        string.ascii_uppercase[13]: ('Percentage sheets:', 0.0, False),
        string.ascii_uppercase[14]: ('Percentage turns:', 0.0, False),
        string.ascii_uppercase[15]: ('Molecular weight:', 0.0, False),
        string.ascii_uppercase[16]: ('Isoelectric point:', 0.0, False),
        string.ascii_uppercase[17]: ('Percentage polar amino acids:', 0.0, False)
    }
    nvdict = vdict.copy()
    # Reference order (alternatively, perhaps easier to keep table order)
    clist = [string.ascii_uppercase[x] for x in [9,8,4,7,12,13,14,15,16,17]]
    # Update score if given and well-formed
    update = False
    if newscore is not None:
        nscore = [x[0] for x in newscore]
        if len( set(nscore) & set(clist) ) > 0:
            for val in vdict:
                nvdict[val] = (vdict[val][0], vdict[val][1], False)
            for val in newscore:
                if val[0] in nvdict:
                    try:
                        newval = float(val[1])
                        nvdict[val[0]] = (nvdict[val[0]][0], newval, True)
                        update = True
                    except:
                        continue
    if not update:
        nvdict = vdict
    score = []
    for x in clist:
        score.append( (x,) + nvdict[x] )
    return score

def updateScore(csvfile, score):
    """ Add or update score column and reorder """
    import string
    head, rows = read_csv(csvfile)
    data = pd_read_csv(csvfile)
    data.index = data.index + 1
    cols = data.columns.tolist()
    sco = pd_Series(np_zeros(len(data[cols[0]])), index=data.index)
    if 'Score' not in cols:
        data['Score'] = sco
        cols = ['Score'] + cols
        data = data[cols]
    colk = list(string.ascii_uppercase)
    for sc in score:
        try:
            coln = colk.index(sc[0])
            val = sc[2]
            checked = sc[3]
            if checked:
                sco += val * data.iloc[:,coln]
        except:
            continue
    data['Score'] = sco
    data = data.sort_values('Score', ascending=False)
    updateMSA(os_path.dirname(csvfile), [[v] for v in data['Seq. ID']])
    data = data.reset_index(drop=True)
    data.index = data.index + 1
    data.rename_axis('Select', axis="columns")
    data.to_csv(csvfile, quoting=csv_QUOTE_ALL, index=False)
    return data

def readRxnCons(consensus):

    f = open(consensus, 'r')

    MnxDir = {}

    for line in f:
        splitdata = line.split()
        Mnx = splitdata[1]
        dirxn = splitdata[2]
        MnxDir[Mnx] = dirxn

    return (MnxDir)

def getMnxSim(rxnInput, datadir, outdir, drxn=0, fp='RDK', pc=None):
    """ Commmand line arguments of quickRsim """

    args = [datadir, fp] + rxnInput + ['-out', os_path.join(outdir,'results_quickRsim.txt')]
    quickRsim.run( quickRsim.arguments(args), pc )
    MnxSim = {}
    MnxDirPref = readRxnCons(os_path.join(datadir, "rxn_consensus_20160612.txt"))
    MnxDirUsed = {}
    EcNumber = {}

    if drxn==1:
        fileout = open(os_path.join(outdir, "results_quickRsim.txt"), 'r')
        for line in fileout:
            splitdata = line.split()
            S1 = splitdata[2]
            S2 = splitdata[3]
            SMILES = splitdata[4]
            EC = ''
            if len(splitdata) > 5:
                EC = splitdata[5]
            Mnx = splitdata[1]
            EcNumber[Mnx] = EC
            try:
                direction = MnxDirPref[Mnx]
                if direction == '1':
                    MnxSim[Mnx]=S1
                    MnxDirUsed[Mnx]= '1'
                elif direction == '-1':
                    MnxSim[Mnx]=S2
                    MnxDirUsed[Mnx]= '-1'
                else:
                    if S1 >= S2:
                        MnxSim[Mnx]=S1
                        MnxDirUsed[Mnx]= '1'
                    elif S2 > S1:
                        MnxSim[Mnx]=S2
                        MnxDirUsed[Mnx]= '-1'
            except KeyError:
                MnxDirPref[Mnx] = "N/A"   # for missing Keys?
                if S1 >= S2:
                    MnxSim[Mnx]=S1
                else:
                    MnxSim[Mnx]=S2
        fileout.close()

        return (MnxSim, MnxDirPref, MnxDirUsed, SMILES, EcNumber)

    else:

        fileout = open(os_path.join(outdir, "results_quickRsim.txt"), 'r')
        for line in fileout:
            splitdata = line.split()
            S1 = splitdata[2]
            S2 = splitdata[3]
            SMILES = splitdata[4]
            EC = ''
            if len(splitdata) > 5:
                EC = splitdata[5]
            Mnx = splitdata[1]
            EcNumber[Mnx] = EC
            if S1 > S2:
                MnxDirUsed[Mnx]='1'
                MnxSim[Mnx]=S1
            elif S2 > S1:
                MnxDirUsed[Mnx]='-1'
                MnxSim[Mnx]=S2
            elif S2 == S1:
                MnxDirUsed[Mnx]='0'
                MnxSim[Mnx]=S2
        fileout.close()

        return (MnxSim, MnxDirPref, MnxDirUsed, SMILES, EcNumber)

def readRxnProp(rxnprop):
    ecrxn = {}
    rxnec = {}
    with open(rxnprop) as handler:
        for line in handler:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            rxnid = row[0]
            ec = row[4].split(';')
            if len(ec) > 0:
                for e in ec:
                    if e == '':
                        continue
                    if e not in ecrxn:
                        ecrxn[e] = set()
                    ecrxn[e].add(rxnid)
                    if rxnid not in rxnec:
                        rxnec[rxnid] = set()
                    rxnec[rxnid].add(e)
    return ecrxn, rxnec

def reactionXref(rxnRefFile, rxnBrenda, rxnSabiork, rxnv3=None):
    rxnref = {}
    for xref in (rxnRefFile, rxnBrenda, rxnSabiork):
        with open(xref) as handler:
            for line in handler:
                if line.startswith('#'):
                    continue
                row = line.rstrip().split('\t')
                rxnref[row[0]] = row[1]
                rxnref['mnx:'+row[1]] = row[1]
    # Downgrade v3 to v2
    if rxnv3 is not None:
        for line in open(rxnv3):
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            try:
                db, rid = row[0].split(':')
                if db == 'deprecated':
                    rxnref['mnx:'+row[1]] = rid
            except:
                continue
    return rxnref

def ecSmiles(ecrxn, rsmi, rxnref):
    smiec = {}
    for ec in ecrxn:
        for r in set(ecrxn[ec]) & set(rsmi):
            s = rsmi[r][0]
            if s not in smiec:
                smiec[s] = set()
            smiec[s].add(ec)
    for ec in ecrxn:
        """ Take the shortest SMARTS being unique for the EC class if possible """
        hits = sorted(set(ecrxn[ec]) & set(rsmi), key=lambda x: len(rsmi[x][0]))
        for h in hits:
            if len(rsmi[h]) > 0:
                rxnid = 'ec:'+ec
                rxnref[rxnid] = h
                if len(smiec[rsmi[h][0]]) == 1:
                       break

def reactionSmiles(rxnSmilesFile):
    rsmi = {}
    with open(rxnSmilesFile) as handler:
        header = handler.readline()
        for line in handler:
            row = line.rstrip().split(',')
            if len(row) > 1:
                rid = row[0]
                rs1 = row[1]
                lr = rs1.split('>>')
                rs2 = lr[1]+'>>'+lr[0]
                rsmi[rid] = (rs1, rs2)

    return rsmi

def pepstats(file, outdir):
    outfile = os_path.join(outdir, "results.pepstats")
    args = ("pepstats -sequence {0} -outfile ".format(file) + outfile)
    os_system(args)

    f = open(outfile, "r")

    hydrop = {}
    weight = {}
    isoelec = {}
    polar = {}

    for line in f:
        if "PEPSTATS of" in line:
            splitdata = line.split()
            seq = splitdata[2].split('_')[0]
        elif "Molecular weight = " in line:
            splitdata = line.split()
            w = splitdata[3]
            weight[seq] = w
        elif "Isoelectric Point = " in line:
            splitdata = line.split()
            i = splitdata[3]
            isoelec[seq] = i
        elif "Polar	" in line:
            splitdata = line.split()
            percent = splitdata[3]
            polar[seq] = percent
            seq
    return (hydrop, weight, isoelec, polar)

def noAmbiguousSeqs(infile, outfile):
    """ Remove ambigous amino acid codes """
    from Bio.Data.IUPACData import protein_letters_1to3, extended_protein_values
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_protein

    newrecords = []
    for record in SeqIO.parse(infile, "fasta"):
        newseq = ''
        for aa in record.seq:
            if aa not in protein_letters_1to3:
                newseq += extended_protein_values[aa][0]
            else:
                newseq += aa
        newrecords.append( SeqRecord(Seq(newseq), id = record.id, description = record.description) )
    SeqIO.write(newrecords, outfile, "fasta")

def garnier(file, outdir):
    fixfile = file+'.fix.fasta'
    noAmbiguousSeqs(file, fixfile)
    outfile = os_path.join(outdir, "garnier.txt")

    args = ("garnier -sequence {0} -outfile ".format(fixfile) + outfile)
    os_system(args)

    f = open(outfile, "r")

    helices = {}
    sheets = {}
    turns = {}
    coils = {}

    for line in f:
        if "Sequence:" in line:
            splitdata = line.split()
            seq = splitdata[2].split('_')[0]
        elif "percent:" in line:
            percents = line.split()
            h = percents[3]
            e = percents[5]
            t = percents[7]
            c = percents[9]
            helices[seq] = h
            sheets[seq] = e
            turns[seq] = t
            coils[seq] = c

    return (helices, sheets, turns, coils)

def updateMSA(outdir, sortrows):
    """ Update the fasta MSA file in the order of the scores """
    align_fasta = os_path.join(outdir, "sequences_aln.fasta")
    if not os_path.exists(align_fasta):
        return
    fasta = {}
    with open(align_fasta) as handler:
        for line in handler:
            if line.startswith('>'):
                splitdata = line.split()
                upid = splitdata[0][1:]
                if len(upid.split('|')) > 1:
                    upid = upid.split('|')[1]
                elif len(upid.split('_')) > 1:
                    upid = upid.split('_')[0]
            if upid not in fasta:
                fasta[upid] = []
            fasta[upid].append(line)

    if len(fasta) == 0:
        return
    outfile = align_fasta
    with open(outfile, 'w') as hw:
        for row in sortrows:
            try:
                pid = row[0]
                if pid in fasta:
                    for seq in fasta[pid]:
                        hw.write(seq)
            except:
                continue

def doMSA(finallistfile, outdir):
    outfile = os_path.join(outdir, "sequences.score_ascii")
    outfile_html = os_path.join(outdir, "sequences.score_ascii.score_html")
    outfile_aln = os_path.join(outdir, "sequences.score_ascii.fasta_aln")
    align_html = os_path.join(outdir, "sequences_score.html")
    align_fasta = os_path.join(outdir, "sequences_aln.fasta")
    treefile = os_path.join(outdir, "sequences.dnd")
    args = ("t_coffee -in {0} -mode quickaln -output=score_ascii,fasta_aln,score_html -outfile ".format(finallistfile) +outfile+ " -newtree "+treefile)
    os_system(args)
    if os_path.exists(outfile_html):
        os_rename(outfile_html, align_html)
    if os_path.exists(outfile_aln):
        os_rename(outfile_aln, align_fasta)

    f = open(outfile, "r")

    cons = {}
    for line in f:
        if "   :  " in line:
            splitdata = line.split()
            upid = splitdata[0]
            if len(upid.split('|')) > 1:
                   upid = upid.split('|')[1]
            elif len(upid.split('_')) > 1:
                   upid = upid.split('_')[0]
            score = splitdata[2]
            cons[upid] = score

    return cons

def read_csv(csvfile):
    rows = []
    if os_path.exists(csvfile):
        with open(csvfile) as handler:
            cv = csv_reader(handler)
            head = next(cv)
            for row in cv:
                rows.append(row)
    return head, rows

def write_csv(csvfilepath, head, rows):
    with open (csvfilepath, 'w') as csvfile:
        writer = csv_writer(csvfile, delimiter=',', quotechar='"', quoting=csv_QUOTE_ALL)
        writer.writerow(head)
        for r in rows:
            writer.writerow(r)

def sort_rows(rows, columns):
    for i in range(len(columns), 0, -1):
        key = columns[i-1]
        if key != 0:
            if key < 0:
                try:
                    rows.sort(key = lambda x: -float(x[abs(key)-1]))
                except:
                    rows.sort(key = lambda x: x[abs(key)-1], reverse=True)
            else:
                try:
                    rows.sort(key = lambda x: float(x[key-1]))
                except:
                    rows.sort(key = lambda x: x[key-1])
    return rows

def write_fasta(fastaFile, targets, pc, short=False, info=False, maxlength=25):
    def shorten(s, maxlength=maxlength):
        ns = re.sub('\s+', '_', s)
        ns = re.sub('\|','_', ns)
        if len(ns) > maxlength:
            ns = ''.join(list(ns)[0:maxlength])+'...'
        return ns
    with open(fastaFile, "w") as f:
        for t in targets:
            try:
                seq = pc.sequence[t]
                if info:
                    desc = ''
                    try:
                        desc = pc.descriptions[t]
                    except:
                        pass
                    org = ''
                    try:
                        org = pc.seqorg[t][1]
                    except:
                        pass
                    ecl = set()
                    if t in pc.UprotToMnx:
                        for r in pc.UprotToMnx[t]:
                            if r in pc.rxnec:
                                ecl |= pc.rxnec[r]
                    ec = ';'.join( sorted(ecl) )
                    st = '_'.join( [t, shorten(desc,30), shorten(org), ec] )
                    print ('>{0} \n{1}'.format(st, seq), file=f)
                elif short:
                    print ('>{0} \n{1}'.format(t, seq), file=f)
                else:
                    fdesc = pc.fulldescriptions[t]
                    print ('>{0} \n{1}'.format(fdesc, seq), file=f)
            except KeyError:
                pass

def short_fasta(fastafile):
    dirname = os_path.dirname(fastafile)
    basename = os_path.basename(fastafile)
    shortfile = 'short_'+basename
    shortname = os_path.join(dirname, shortfile)
    with open(fastafile) as fasta, open(shortname, 'w') as fastaw:
        for line in fasta:
            if line.startswith('>'):
                head = line.split('OS=')
                if len(head[0].split('|')) >= 2:
                    seqid = head[0].split('|')[1]
                    line = '>'+seqid+'\n'
            fastaw.write(line)
    return shortfile

def extend_sequences(initialfastafile, fastafile, workfolder, noMSA):
    """ Extend the fasta file """
    csvfile = os_path.join(workfolder, 'selenzy_results.csv')
    shortfile = short_fasta(os_path.join(workfolder, fastafile))
    try:
        # TO do: Check that does not assume Uniprot format...
        (sequence, names, descriptions,
         fulldescriptions, osource, pexistence) = readFasta(workfolder, fastafile, limit=1000)

        if len(sequence) == 0:
            return csvfile

        (hydrop, weight, isoelec, polar, helices,
         sheets, turns, coils) = sequence_properties(os_path.join(workfolder, shortfile))
        fasta = open(os_path.join(workfolder, initialfastafile)).readlines()
        fasta.extend(open(os_path.join(workfolder, fastafile)).readlines())
        with open(os_path.join(workfolder, initialfastafile), 'w') as handler:
            for line in fasta:
                handler.write(line)
        shortfile = short_fasta(os_path.join(workfolder, initialfastafile))
        if not noMSA:
            cons = conservation_properties(os_path.join(workfolder, shortfile))
        else:
            cons = ({},{})

    except:
        return csvfile
    # Extend the csvfile
    head, rows = read_csv(csvfile)
    # Update the conservation scores
    if not noMSA:
        for r in rows:
            seqid = r[head.index('Seq. ID')]
            if seqid in cons:
                r[head.index('Consv. Score')] = cons[seqid]
    for k in range(0, len(fulldescriptions)):
        n = names[k]
        try:
            conservation = cons[n]
        except:
            conservation = 0
        try:
            h = helices[n]
            e = sheets[n]
            t = turns[n]
            c = coils[n]
        except:
            h = '-'
            e = '-'
            t = '-'
            c = '-'
        try:
            w = weight[n]
            i = isoelec[n]
            pol = polar[n]
        except:
            w = '-'
            i = '-'
            pol = '-'
        try:
            description = descriptions[n]
            source = osource[n]
            ext = pexistence[n]
        except:
            description = '-'
            source = '-'
        row = [0.0, n, description, source, -1,
               '-', '-', ext, conservation, 1, 1, 1,
               h, e, t, c, w, i, pol, '-', '-']
        rows.append( row )
    write_csv(csvfile, head, rows)
    return csvfile

def sequence_properties(fastaShortNameFile):
    #analysis of final list of sequences
    (hydrop, weight, isoelec, polar) = pepstats(fastaShortNameFile, os_path.dirname(fastaShortNameFile))
    (helices, sheets, turns, coils) = garnier(fastaShortNameFile,  os_path.dirname(fastaShortNameFile))
    return hydrop, weight, isoelec, polar, helices, sheets, turns, coils

def conservation_properties(fastaFile):
    cons = doMSA(fastaFile,  os_path.dirname(fastaFile))
    return cons

def analyse(rxnInput, targ, datadir, outdir, csvfilename, pdir=0, host='83333', NoMSA=False, pc=None, fp='RDK'):


    datadir = os_path.join(datadir)
    outdir = os_path.join(outdir)

#    rxnname = os_path.splitext(rxn)[0]
#    csvname = rxnname.rsplit('/', 1)[-1]

    if csvfilename:
        csvfilename = csvfilename
    else:
        csvfilename = "results_selenzy.csv"
    if not pc:
        pc = preLoad(datadir)

    print ("Running quickRsim...", end='')
    try:
        (MnxSim, MnxDirPref, MnxDirUsed, Smiles, EcNumber) = getMnxSim(rxnInput, datadir, outdir, pdir, fp, pc)
        print_OK()
    except:
        print_FAILED()
        return False, pc


    sequence = pc.sequence
    names = pc.names
    descriptions = pc.descriptions
    fulldescriptions = pc.fulldescriptions
    osource = pc.osource
    pexistence = pc.pexistence
    seqorg = pc.seqorg
    tax = pc.tax
    MnxToUprot = pc.MnxToUprot
    upclst = pc.upclst
    clstrep = pc.clstrep
    smir = pc.smir

    list_mnx = sorted(MnxSim, key=MnxSim.__getitem__, reverse=True)  #allow user to manipulate window of initial rxn id list
    print ("Creating initial MNX list...")
    targets = set()
    UprotToMnx = {}


    # first creating fasta file, f, for further data extraction
    for x in list_mnx:
        up = MnxToUprot.get(x)
        if up is not None:
            for y in up:
                if len(targets) >= int(targ):    # allow user to manipulate desired number of entries for resulting table
                    break
                else:
                    UprotToMnx[y] = x
                    try:
                        targets.add(y)
                    except KeyError:
                        pass
            else:
                continue
            break
    fastaFile = os_path.join(outdir, "sequences.fasta")
    write_fasta(fastaFile, targets, pc)
    # Avoid issues with sequence ids
    fastaShortNameFile = os_path.join(outdir, "seqids.fasta")
    fastaInfoNameFile = os_path.join(outdir, "seqinfo.fasta")
    write_fasta(fastaShortNameFile, targets, pc, short=True)
    write_fasta(fastaInfoNameFile, targets, pc, info=True)


    (hydrop, weight, isoelec, polar, helices, sheets, turns, coils) = sequence_properties(fastaShortNameFile)
    if not NoMSA:
        cons = conservation_properties(fastaInfoNameFile)
    else:
        cons = {}

    print ("Acquiring sequence properties...")
    # final table, do all data and value storing before this!
    tdist = {}
    rows = []


    for y in targets:
        try:
            # Essential sequence information
            desc = descriptions[y]
            fdesc = fulldescriptions[y]
            org = osource[y]
            ext = pexistence[y]
            mnx = UprotToMnx[y]
            rxnsimpre = float(MnxSim[mnx])
            rxnsim = float("{0:.5f}".format(rxnsimpre))

            # Non-essential sequence information
            try:
                cn = upclst.get(y)
                repid = clstrep[cn]
            except:
                cn = 0
                repid = y
            try:
                conservation = float(cons[y])
            except:
                conservation = 0.0
            try:
                ecid = EcNumber[mnx]
            except:
                ecid = ''
            try:
                h = helices[y]
                e = sheets[y]
                t = turns[y]
                c = coils[y]
            except:
                h = '-'
                e = '-'
                t = '-'
                c = '-'
            try:
                w = weight[y]
                i = isoelec[y]
                pol = polar[y]
            except:
                w = '-'
                i = '-'
                pol = '-'
            try:
                rxndirused = MnxDirUsed[mnx]
                rxndirpref = MnxDirPref[mnx]
            except:
                rxndirused = 1
                rxndirpred = 1

            mnxSmiles = ''
            if mnx in smir:
                if rxndirused == 1:
                    mnxSmiles = smir[mnx][0]
                else:
                    mnxSmiles = smir[mnx][1]
            if org not in tdist:
                if y in seqorg:
                    tdist[org] = taxDistance(tax, host, seqorg[y][0])
                else:
                    tdist[org] = '-'

            rows.append( (y, desc, org, tdist[org], mnx, ecid, ext, conservation, rxnsim, rxndirused, rxndirpref, h, e, t, c, w, i, pol, Smiles, mnxSmiles) )


        except KeyError:
            pass
    sortrows = sort_rows(rows, (-10, -9, 4) )
    updateMSA(outdir, sortrows)

    head = ('Seq. ID','Description', 'Organism Source', 'Tax. distance', 'Rxn. ID', 'EC Number', 'Uniprot protein evidence', 'Consv. Score',
            'Rxn Sim.', "Direction Used", "Direction Preferred",
            '% helices', '% sheets', '% turns', '% coils', 'Mol. Weight', 'Isoelec. Point', 'Polar %','Query', 'Hit')

    write_csv(os_path.join(outdir, csvfilename), head, sortrows)

    print ("CSV file created.")
    return True, pc

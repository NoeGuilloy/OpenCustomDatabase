from OpenCMS.utils import *
from OpenCMS.SNP_handle import *
OP_protein_fasta = '/home/noeguill/OpenCMS/Databases/OpenProt_DB_Ensemble_1_6.fasta'
transcrit_fasta = '/home/noeguill/OpenCMS/Databases/gencode.v29.transcripts.fa'
OP_tsv = '/home/noeguill/OpenCMS/Databases/human-openprot-r1_6-refprots+altprots+isoforms-ensembleonly.tsv'

class OpenCMS:
    def __init__(self, 
    vcf_path,
    expname,
    input_kallisto = None,
    trxnumber = None,          #10000-1000000
    tpmnumber = None,          #0-1000000
    ipban = None,              #yes
    trxsave = None,
    trxexclude = None,
    annotation='ensembl'):
        self.vcf_path = vcf_path
        self.expname = expname
        self.input_kallisto = input_kallisto
        self.trxnumber = trxnumber
        self.tpmnumber = tpmnumber
        self.annotation = annotation
        self.ipban = ipban
        self.trxsave = trxsave
        self.trxexclude = trxexclude

    def run(self, verbose=True):
        if verbose:
            print('checking files')
        checkex = checking_trx_files(self.trxexclude)
        checksv = checking_trx_files(self.trxsave)
        checkka = abundance_check(self.input_kallisto)
        if checkex == False:
            return print('Make sure your exclusion transcrit file is well written')
        if checksv == False:
            return print('Make sure your save transcrit file is well written')
        if checkka == False:
            return print('Make sure your kallisto-quant file is not corrupted')
        print('running Openvar...')
        parsed_snpeff = OpenVar_analysis(self.vcf_path, self.expname)
        print('Parsing...')
        protvariantfile = get_protvcf_file(parsed_snpeff, self.expname, self.vcf_path)
        print('..Done')
        print('calculs...')
        prot_syno = get_synonyms_prot(OP_protein_fasta)
        start_codon = get_start_codon(OP_tsv, transcrit_fasta)
        fasta_dict = get_fasta_dict(OP_protein_fasta)
        var_by_prot,transcrit_prot = parse_protvcf_file(protvariantfile)
        seqname_seq = get_all_mut_sequences(var_by_prot,transcrit_prot,start_codon,fasta_dict,prot_syno,self.ipban)
        print('phase2')
        if self.input_kallisto:
            print('started')
            trx_allprot= append_wt_prot_to_transcrit_by_fasta(OP_protein_fasta,self.ipban)
            print('append_wt_prot_to_transcrit_ENS done')
            trx_allprot= get_mutated_protbytranscrit(seqname_seq,transcrit_prot,trx_allprot)
            print('get_mutated_protbytranscrit done')
            AllProtInMyDB,  effective_threshold = get_100_prot(self.input_kallisto, prot_syno,trx_allprot, self.trxnumber, self.trxsave, self.tpmnumber, self.trxexclude)
            DB_custom = assembling_headers_sequences(AllProtInMyDB,seqname_seq,prot_syno,fasta_dict)
            DB_custom = remove_duplicata_from_db(DB_custom)
            write_Fasta_DB(DB_custom, self.expname, self.vcf_path,effective_threshold)

def abundance_check(input_kallisto):
    if input_kallisto:
        with open(input_kallisto, 'r') as f:
            for n,l in enumerate(f):
                if n == 0:
                    if l != 'target_id\tlength\teff_length\test_counts\ttpm\n':
                        return False
                else:
                    if len(l.split('\t'))!=5:
                        return False
    return True

def checking_trx_files(trxfile):
    if trxfile:
        with open(trxfile, 'r') as f:
            for n,l in enumerate(f):
                if l[0:4] != 'ENST' and  l[0:3] != 'NM_' and  l[0:3] != 'XM_':
                    return False
                if len(l)>40:
                    return False
    return True

def get_all_mut_sequences(var_by_prot,transcrit_prot,start_codon,fasta_dict,prot_syno,ipban):
    seqname_seq=dict()
    for acc,svar in var_by_prot.items():
        if ipban=='yes':
            if acc[0:3]=='II_' or acc [0:3]=='IP_':continue
        regroupement_HGVS_C = list()
        regroupement_HGVS_P = list()
        for var in svar:
            parsed_HGVS_C = Parse_HGVS_C(var['HGVS_C'])
            parsed_HGVS_P = var['HGVS_P'].split('.',1)[1]
            regroupement_HGVS_C.extend(parsed_HGVS_C)
            regroupement_HGVS_P.append(parsed_HGVS_P)
        sorted_HGVS_C = sort_sequences(regroupement_HGVS_C)
        seqname = acc+'@'+''.join(regroupement_HGVS_P)
        mutated_sequence = modify_transcript_sequence(sorted_HGVS_C,acc,transcrit_prot,start_codon)
        translated_mutated_sequence = translate(mutated_sequence)
        if 'synonymous_variant' in mutated_sequence or translated_mutated_sequence == 'start_lost' or len(translated_mutated_sequence)<7:continue
        elif translated_mutated_sequence in seqname_seq.values():continue
        else:
            seqname_seq[seqname]=translated_mutated_sequence
    seqname_seq=remove_fakevariant(seqname_seq,fasta_dict,prot_syno)
    return seqname_seq

def append_wt_prot_to_transcrit_by_fasta (fasta,ipban):
    trx_allprot = defaultdict(list)
    with open(fasta, 'r')as f:
        for n,l in enumerate(f):
            if '>' not in l:continue
            parsedheader=parse_fasta_header(l)
            trxs=parsedheader['TA'].split(',')
            prxs=parsedheader['PA'].split(',')
            for trx in trxs:
                for prx in prxs:
                    if ipban == "yes":
                        if prx[0:3]=='II_' or prx [0:3]=='IP_':continue
                    trx_allprot[trx.split('.')[0]].append(prx.split('.')[0])
    return trx_allprot

def get_100_prot(trx_expression,prot_syno,trx_allprot,trxnumber,trxsave,tpmnumber,trxexclude):
    trx_expression_sorted = get_trxs_by_tpm_from_kallisto(trx_expression,trxexclude)
    AllProtInMyDB = set()
    
    if trxnumber:
        treshold = trxnumber
    else:
        treshold = 100000
    for prot,tpm in trx_expression_sorted.items():
        if tpmnumber:
            if tpm>tpmnumber:
                for y in trx_allprot[prot]:    
                    if y in prot_syno:
                        if prot_syno[y] not in AllProtInMyDB:
                            AllProtInMyDB.add(prot_syno[y])
                    else:
                        AllProtInMyDB.add(y)
                    if y not in AllProtInMyDB and prot_syno[y] not in AllProtInMyDB:
                        print('not_taken',prot,y)
            else:
                print(prot,tpm)
                return(AllProtInMyDB,tpm)
        else:
            for y in trx_allprot[prot]:    
                if len(AllProtInMyDB)<treshold:
                    if y in prot_syno:
                        if prot_syno[y] not in AllProtInMyDB:
                            AllProtInMyDB.add(prot_syno[y])
                    else:
                        AllProtInMyDB.add(y)
                    if y not in AllProtInMyDB and prot_syno[y] not in AllProtInMyDB:
                        print('not_taken',prot,y)
                else:
                    print(prot,tpm)
                    return(AllProtInMyDB,tpm)
    if trxsave:
        with open(trxsave, 'r') as f:
            for n,l in enumerate(f):
                for y in trx_allprot[l.split('.')[0]]:
                    if y in prot_syno:
                        if prot_syno[y] not in AllProtInMyDB:
                            AllProtInMyDB.add(prot_syno[y])
                    else:
                        AllProtInMyDB.add(y)
    print(prot,tpm)            
    return(AllProtInMyDB,tpm)

def get_trxs_by_tpm_from_kallisto (trx_expression,trxexclude):
    if trxexclude:
        liste_exclusion = list()
        with open(trxexclude, 'r') as fex:
            for n,l in enumerate(fex):
                liste_exclusion.append(l.split('\n')[0].split('.')[0])
    else :
        liste_exclusion = " "
        
    trxp_tpms = {}
    with open(trx_expression, 'r') as f:
        for n,l in enumerate(f):
            ls = l.strip().split('\t')
            if n==0:
                keys = ls
                continue
            line = dict(zip(keys, ls))
            trxp = line['target_id'].split('.')[0]
            tpm = float(line['tpm'])
            if tpm >0:
                if trxp in liste_exclusion:continue
                trxp_tpms[trxp] = tpm
    trxp_tpms_sorted_pour_tresh = OrderedDict(sorted(trxp_tpms.items(), key=itemgetter(1), reverse=True))
    return (trxp_tpms_sorted_pour_tresh)

def get_mutated_protbytranscrit(seqname_seq,transcrit_prot,trx_allprot):
    for prot_mut in seqname_seq:
        for prot,transcrit in transcrit_prot.items():
            if prot in prot_mut:
                trx_allprot[transcrit].append(prot_mut)
    return trx_allprot

def modify_transcript_sequence(HGVS_C_snp_sorted,protacc,transcrit_prot,start_codon):
    trx = str(transcrit_prot[protacc])
    acc = protacc.split('.',1)[0]
    mseq = [x for x in str(start_codon[acc+'^'+trx])]
    for variant in HGVS_C_snp_sorted:
        if variant['effect'] == 'del':
            if mseq[int(variant['pos'])-1]:
                mseq[int(variant['pos'])-1] = ""
        if variant['effect'] == 'ins':
            mseq.insert(int(variant['pos']), variant['nt'])
        if variant['effect'] == 'replace':
            mseq[int(variant['pos'])-1] = variant['nt']
    mseq[:] = [x for x in mseq if x !='']
    if translate(''.join(mseq)) == translate(''.join([x for x in str(start_codon[acc+'^'+trx])])):
        return 'synonymous_variant',acc
    if len(translate(''.join(mseq))) < 8:
        return 'toosmall',acc
    mseq = "".join(mseq)
    return mseq

def get_m_wtprot(seqname_seq):
    m_wtprot = list()
    for prot_mut in seqname_seq:
        m_wtprot.append(prot_mut.split('@')[0])
    return m_wtprot

def assembling_headers_sequences(AllProtInMyDB,Msequence,prot_syno,fasta_dict):
    DB_custom = dict()
    for acc in AllProtInMyDB:
        regular_acc = acc.split('@')[0].split('.')[0]

        if acc in Msequence:
            if regular_acc in prot_syno:
                header = acc+'|'+str(fasta_dict[prot_syno[regular_acc]].description).split('|')[1]
                DB_custom[header] = Msequence[acc]
            else:
                header = acc+'|'+str(fasta_dict[regular_acc].description).split('|')[1]
                DB_custom[header] = Msequence[acc]
        else:
            if acc in prot_syno:
                header = acc+'|'+str(fasta_dict[prot_syno[regular_acc]].description).split('|')[1]
                DB_custom[header] = str(fasta_dict[prot_syno[regular_acc]].seq)
            else:
                header = acc+'|'+str(fasta_dict[regular_acc].description).split('|')[1]
                DB_custom[header] = str(fasta_dict[regular_acc].seq)
    
    return DB_custom

def stat_summary(effective_threshold,DB_custom,expname,vcf_path):
    Number_IP = 0
    Number_IPvar = 0
    Number_II = 0
    Number_IIvar = 0
    Number_ref = 0
    Number_refvar = 0
    for acc in DB_custom.keys():
        if acc[:3]=='IP_':
            if '@' in acc:
                Number_IPvar = Number_IPvar+1
            else:
                Number_IP = Number_IP+1
        elif acc[:3]=='II_':
            if '@' in acc:
                Number_IIvar = Number_IIvar+1
            else:
                Number_II = Number_II+1      
        else:
            if '@' in acc:
                Number_refvar = Number_refvar+1
            else:
                Number_ref = Number_ref+1
    filename = vcf_path.split('/')[-1]
    path = filename.replace('.vcf','')+'_result/'+expname+'summary.tsv'
    with open(path, 'w') as f:
        f.write('effective_threshold\taltProt (IP) \tnovel isoform (II)\trefprot\taltProt_variants (IP) \tnovel isoform_variants (II)\trefprot_variants\n')
        f.write(str(effective_threshold)+'\t'+str(Number_IP)+'\t'+str(Number_II)+'\t'+str(Number_ref)+'\t'+str(Number_IPvar)+'\t'+str(Number_IIvar)+'\t'+str(Number_refvar)+'\n')

def write_Fasta_DB(DB_custom,expname,vcf_path,effective_threshold):
    filename = vcf_path.split('/')[-1]
    path = filename.replace('.vcf','')+'_result/'+expname+'.fasta'
    stat_summary(effective_threshold,DB_custom,expname,vcf_path)
    with open(path, 'w') as f:
        for acc, prot_seq in DB_custom.items():
            f.write('>'+acc+'\n')
            f.write(truncate(prot_seq+'\n'))
    return print(path+'is done')

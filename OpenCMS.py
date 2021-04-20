from OpenCMS.utils import *
from OpenCMS.SNP_handle import *
OP_protein_fasta = './Databases/OpenProt_DB_Ensemble_1_6.fasta'
transcrit_fasta = './Databases/gencode.v29.transcripts.fa'
OP_tsv = './Databases/human-openprot-r1_6-refprots+altprots+isoforms-+uniprot2019_03_01.tsv'

class OpenCMS:
    def __init__(self, 
    vcf_path, 
    output_db_path,
    input_kallisto,
    expname, 
    annotation='ensembl'):
        self.vcf_path = vcf_path
        self.expname = expname
        self.output_db_path = output_db_path
        self.input_kallisto = input_kallisto
        self.annotation=annotation

    def run(self, verbose=True):
        if verbose:
            print('running Openvar...')
        parsed_snpeff = OpenVar_analysis(self.vcf_path, self.expname)
        print('Parsing...')
        protvariantfile = get_protvcf_file(parsed_snpeff,expname)
        print('..Done')
        print('calculs...')
        prot_syno = get_synonyms_prot(OP_protein_fasta)
        start_codon = get_start_codon(OP_tsv, transcrit_fasta)
        fasta_dict = get_fasta_dict(OP_protein_fasta)
        var_by_prot,transcrit_prot = parse_protvcf_file(protvariantfile)
        seqname_seq = self.get_all_mut_sequences(var_by_prot,transcrit_prot,start_codon,fasta_dict,prot_syno,ipban='')


    def get_all_mut_sequences(self,var_by_prot,transcrit_prot,start_codon,fasta_dict,prot_syno,ipban):
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


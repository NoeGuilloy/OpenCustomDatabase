from utils import *
from SNP_handle import *


class OpenCMS:
    def __init__(self, 
    input_vcf, 
    output_db_path,
    input_kallisto, 
    annotation='ensembl'):
        self.input_vcf = input_vcf
        self.output_db_path = output_db_path
        self.input_kallisto = input_kallisto
        self.annotation=annotation

    def run(self, verbose=True):
        if verbose:
            print('running snpeff...')
        parsed_snpeff = self.snpeff()
        
        modifed_tab = utils.compute_mod()
        self.get_mutated_seq(parsed_snpeff)

        prot_syno = OpenCMS_utils.get_synonyms_prot(protein_fasta)
        start_codon = OpenCMS_utils.get_start_codon(tsv, transcrit_fasta)
        fasta_dict = OpenCMS_utils.get_fasta_dict(protein_fasta)
        protvariantfile = OpenCMS_utils.get_protvcf_file(parsenpff,expname)
        var_by_prot,transcrit_prot = OpenCMS_utils.parse_protvcf_file(protvariantfile)
        seqname_seq = get_all_mut_sequences(var_by_prot,transcrit_prot,start_codon,fasta_dict,prot_syno,ipban)


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

def get_mutated_seq(self, tabfile):
        return
        
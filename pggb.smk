

def get_hg_filter_ani(wildcards):
    if wildcards.p >= 90:
        return 30
    elif wildcards.p >= 80:
        return 10
    else:
        return 0

rule wfmash_mapping:
    input:
        fasta = ''
    output:
        mapping = 'pggb/p{p}_s{segment_length}/mapping.paf'
    params:
        hg_ani_diff = get_hg_filter_ani,
        block_length = lambda wildcards: wildcards.segment_length * 3
    shell:
        '''
        wfmash \
          -s {wildcards.segment_length} \
          -l {params.block_length} \
          -p {wildcards.p} \
          -n 1 \
          -k 19 \
          -H 0.001 \
          -Y "#" \
          -t {threads} \
          --tmp-base $TMPDIR \
          --hg-filter-ani-diff {params.hg_ani_diff} \
          {input.fasta} \
          --lower-triangular \
          --approx-map \
          > {output.mapping}
         '''

rule wfmash_alignment:
    input:
        fasta = '',
        mapping = rules.wfmash_mapping.output['mapping']
    output:
        'pggb/p{p}_s{segment_length}/alignment.paf
    params:
        hg_ani_diff = get_hg_filter_ani,
        block_length = lambda wildcards: wildcards.segment_length * 3
    shell:
        '''
        wfmash \
          -s {wildcards.segment_length} \
          -l {params.block_length} \
          -p {wildcards.p} \
          -n 1 \
          -k 19 \
          -H 0.001 \
          -Y "#" \
          -t {threads} \
          --tmp-base $TMPDIR \
          {input.fasta} \
          --lower-triangular \
          --hg-filter-ani-diff {params.hg_ani_diff}\
          -i {input.mapping} \
        > {output.alignment}
        '''

rule seqwish:
    input:
        fasta = '',
        alignment = rules.wfmash_alignment.output['paf']
    output:
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.seqwish.gfa
    shell:
        '''
            seqwish \
                -s {input.fasta} \
                -p {input.alignment} \
                -k {wildcards.k} \
                -f 0 \
                -g {output.gfa} \
                -B 10000000 \
                -t {threads} \
                --temp-dir $TMPDIR \
                -P \
        '''

def POA_params(wildcards.?):
    match wildcards.POA_setting:
        case 'asm5':
            return "1,19,39,3,81,1"
        case 'asm10':
            return "1,9,16,2,41,1"
        case 'asm15':
            return "1,7,11,2,33,1"
        case 'asm20':
            return "1,4,6,2,26,1"

rule smoothxg:
    input:
        fasta = '',
        gfa = rules.seqwish.output['gfa']
    output:
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.seqwish.gfa
    params:
        block_id_min = lambda wildcards: round(wildcards.p / 4,4),
        n_haps = lambda wildcards, input: sum(1 for _ in open(input.fasta[2])),
        POA_pad_depth = '"$pad_max_depth * $n_haps" | bc)'
        POA_lengths = '700,900,1100',
        POA_params = POA_params()
    shell:
        '''
            smoothxg \
             -t {threads} \
             -T {threads} \
             -g {input.gfa} \
             -r {params.n_haps} \
             --base $TMPDIR \
             --chop-to 100 \
             -I {params.block_id_min} \
             -R 0 \
             -j 0 \
             -e 0 \
             -l {params.POA_lengths} \
             -P {params.POA_params} \
             -O 0.001 \
             -Y {params.POA_pad_depth} \
             -d 0 -D 0 \
             -Q "$consensus_prefix" \
             $consensus_params \
             -o {output.gfa}
        '''

rule gffafix:
    input:
        gfa = rules.smoothxg.output['gfa']
    output:
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.gffafix.gfa'
        affixes = 'pggb/p{p}_s{segment_length}/k{k}.gffafix.affixes.tsv.gz'
    shell:
        '''
        gfaffix {input.gfa} -o {output.gfa} |\
        pigz -p 2 -c > {output.affixes}
        '''

#Is this helpful?
rule vg_path_normalise:
    input:
        ''
    output:
        ''
    shell:
        '''
        vg mod -X 1024 -
        vg paths -x - -n -Q options.reference[0] -t {threads}
        '''

rule odgi_unchop:
    input:
        gfa = rules.gffafix.output['gfa']
    output:
        og = 'pggb/p{p}_s{segment_length}/k{k}.unchop.og',
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.unchop.gfa',
    shell:
        '''
        odgi build -t {threads} -P -g {input.gfa} -o - -O |\
        odgi unchop -P -t {threads} -i - -o - |\
        odgi sort -P -p Ygs --temp-dir $TMPDIR -t {threads} -i - -o - |\
        tee {output.og} |\
        odgi view -i - -g > {output.gfa}
        '''

rule all:
    input:
        'pggb/p80_s5000/k19.POAasm5.unchop.og'

def get_hg_filter_ani(wildcards):
    if float(wildcards.p) >= 90:
        return 30
    elif float(wildcards.p) >= 80:
        return 10
    else:
        return 0

ruleorder: split_approx_mappings_in_chunks > wfmash

wildcard_constraints:
    mode = r'mapping|alignment',
    chunk = r'\.\d+|'

rule wfmash:
    input:
        fasta = 'test.fa',
        mapping = lambda wildcards: 'pggb/p{p}_s{segment_length}/mapping{chunk}.paf' if wildcards.mode == 'alignment' else []
    output:
        paf = 'pggb/p{p}_s{segment_length}/{mode}{chunk,(\.\d+|)}.paf'
    params:
        hg_ani_diff = get_hg_filter_ani,
        block_length = lambda wildcards: wildcards.segment_length * 3,
        mapping = lambda wildcards, input: f'-i {input.mapping}' if wildcards.mode == 'alignment' else ''
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
          {params.mapping} > {output.paf}
         '''

#based off https://github.com/waveygang/wfmash/blob/main/scripts/split_approx_mappings_in_chunks.py
rule split_approx_mappings_in_chunks:
    input:
        mapping = expand(rules.wfmash.output['paf'],mode='mapping',chunk='',allow_missing=True)
    output:
        paf = expand('pggb/p{p}_s{segment_length}/mapping.{chunk}.paf',chunk=range(config.get('wfmash_chunks',1)),allow_missing=True)
    run:
        rank_to_mapping_dict, mapping_list = {}, []

        with open(input['mapping']) as f:
            for rank, line in enumerate(f):
                # We could avoid keeping everything in memory by reading the file again later
                rank_to_mapping_dict[rank] = line

                _, _, query_start, query_end, _, _, _, target_start, target_end, _, _, _, estimated_identity = line.strip().split('\t')[:13]

                num_mapped_bases = max(int(query_end) - int(query_start), int(target_end) - int(target_start))
                estimated_identity = float(estimated_identity.split('id:f:')[1]) / 100.0

                # High divergence makes alignment more difficult
                weight = num_mapped_bases * (1 - estimated_identity)

                mapping_list.append((rank, weight))

        chunk_list = [[] for i in range(num_of_chunks)]
        sums = [0] * num_of_chunks
        i = 0
        for e in mapping_list:
            chunk_list[i].append(e)
            sums[i] += e[1]
            i = sums.index(min(sums))

        # Collect the ranks from the tuples to generate balanced chunks
        for num_chunk, element_list in enumerate(chunk_list):
            with open( f'.chunk_{num_chunk}.paf', 'w') as fw:
                for rank, _ in element_list:
                    fw.write(rank_to_mapping_dict[rank])


rule wfmash_concat:
    input:
        expand(rules.wfmash.output['paf'],mode='alignment.',chunk=range(config.get('wfmash_chunks',1)),allow_missing=True)
    output:
        paf = 'pggb/p{p}_s{segment_length}/alignment.concat.paf'
    localrule: True
    shell:
        '''
        cat {input} > {output}
        '''

rule seqwish:
    input:
        fasta = 'test.fa',
        alignment = expand(rules.wfmash_concat.output['paf'],allow_missing=True) if config.get('wfmash_chunks',1) > 1 else expand(rules.wfmash.output['paf'],mode='alignment',chunk='',allow_missing=True)
    output:
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.seqwish.gfa'
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
                -P
        '''

def POA_params(wildcards):
    match wildcards.POA:
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
        fasta = multiext('test.fa.gz','.fai','.gzi'),
        gfa = rules.seqwish.output['gfa']
    output:
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.POA{POA}.smoothxg.gfa'
    params:
        block_id_min = lambda wildcards: round(float(wildcards.p) / 4,4),
        n_haps = lambda wildcards, input: sum(1 for _ in open(input.fasta[1])),
        POA_pad_depth = lambda wildcards, input: 100 * sum(1 for _ in open(input.fasta[1])),
        POA_lengths = '700,900,1100',
        POA_params = POA_params
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
             -o {output.gfa}
        '''

rule gffafix:
    input:
        gfa = rules.smoothxg.output['gfa']
    output:
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.POA{POA}.gffafix.gfa',
        affixes = 'pggb/p{p}_s{segment_length}/k{k}.POA{POA}.gffafix.affixes.tsv.gz'
    shell:
        '''
        gfaffix {input.gfa} -o {output.gfa} |\
        pigz -p 2 -c > {output.affixes}
        '''

#Is this helpful?
rule vg_path_normalise:
    input:
        fasta = multiext('test.fa.gz','.fai','.gzi'),
        gfa = rules.gffafix.output['gfa']
    output:
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.POA{POA}.vg.gfa'
    params:
        reference = lambda wildcards, input: open(input.fasta[1]).readline().rstrip()
    shell:
        '''
        vg convert -g {input.gfa} -t {threads} -x |\
        #vg mod -X 1024 {input.gfa} |\ #unneeded if already unchopped
        vg paths -x - -n -Q {params.reference} -t {threads}
        '''

rule odgi_unchop:
    input:
        gfa = rules.vg_path_normalise.output['gfa']
    output:
        og = 'pggb/p{p}_s{segment_length}/k{k}.POA{POA}.unchop.og',
        gfa = 'pggb/p{p}_s{segment_length}/k{k}.POA{POA}.unchop.gfa',
    shell:
        '''
        odgi build -t {threads} -P -g {input.gfa} -o - -O |\
        odgi unchop -P -t {threads} -i - -o - |\
        odgi sort -P -p Ygs --temp-dir $TMPDIR -t {threads} -i - -o - |\
        tee {output.og} |\
        odgi view -i - -g > {output.gfa}
        '''

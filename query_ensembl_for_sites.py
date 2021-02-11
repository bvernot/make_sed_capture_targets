import os,sys
from ensembldb3 import HostAccount, Compara

# compara = Compara(['human', 'chimp', 'macaque'], release=75)

offset = int(sys.argv[1])
offset = offset + 1 # account for 1bp difference btwn annotations?  this is needed to make a 25bp offset w/ 52bp len produce 25X26
length = int(sys.argv[2])


alleles_file = open(sys.argv[3])
# [bionc02 soil_capture_2017/site_categories_for_capture]$ head admixture_array/admixture_array.snp.alleles
# 1       847983  C       T
# 1       853089  G       C
# 1       853596  A       G
# 1       854793  A       C
# 1       867151  G       A

allelic_state = {}
for line in alleles_file:
    chrom, pos, a1, a2 = line.rstrip().split()
    pos = int(pos)
    allelic_state[(chrom,pos)] = (a1,a2)
    pass

print('.. read %d allelic states' % len(allelic_state))


map_rpt_file = open(sys.argv[4])
# head reich/reich_positions.snp.ape.bi.tv.rpt_map
# 1       842012  842013  0       100     1
# 1       846863  846864  0       104     1
# 1       869302  869303  10      88      1

map_states = {}
for line in map_rpt_file:
    chrom, _, pos, repeat_mask_bases, heng99_bases, heng99_target_overlap = line.rstrip().split()
    pos = int(pos)
    repeat_mask_bases = int(repeat_mask_bases)
    heng99_bases = int(heng99_bases)
    heng99_target_overlap = int(heng99_target_overlap)
    map_states[(chrom,pos)] = (repeat_mask_bases, heng99_bases, heng99_target_overlap)
    pass

print('.. read %d map/repeat states' % len(map_states))


ape_file = open(sys.argv[5])
# head reich/reich_positions.snp.ape
# 1       842013  T       T       T       T       N       T       T       T       T
# 1       846864  G       G       G       G       G       G       G       G       C
# 1       869303  C       C       C       C       C       C       C       C       C

## columns:
# 1 chr
# 2 pos
# 3 hg19
# 4 pantro4
# 5 panpan1.1
# 6 gorgor3
# 7 ponabe2
# 8 rhemac3
# 9 Vindija33.19 (with ambiguity codes for hets)
# 10 Altai (with ambiguity codes for hets)
# 11 Denisova (with ambiguity codes for hets)

ape_states = {}
for line in ape_file:
    chrom, pos, hg19, pantro, panpan, gorgor, ponabe, rhemac, vind, altai, den = line.rstrip().split()[:11]
    pos = int(pos)
    ape_states[(chrom,pos)] = (hg19, pantro, panpan, gorgor, ponabe, rhemac, vind, altai, den)
    pass

print('.. read %d ape states' % len(ape_states))

 


spc_all = ['callithrix_jacchus',
           'bos_taurus',
           'pan_troglodytes',
           'equus_caballus',
           'ovis_aries',
           'canis_familiaris',
           'macaca_mulatta',
           'mus_musculus',
           'oryctolagus_cuniculus',
           'gorilla_gorilla',
           'sus_scrofa',
           'felis_catus',
           'pongo_abelii',
           'rattus_norvegicus',
           'homo_sapiens']

spc_chimp = set(('pan_troglodytes',))
spc_apes = set(('pongo_abelii', 'gorilla_gorilla', 'pan_troglodytes'))
spc_monkey = set(('macaca_mulatta', 'callithrix_jacchus', 'chlorocebus_sabaeus', 'papio_anubis'))
spc_rodents = set(('rattus_norvegicus', 'mus_spretus_spreteij', 'mus_musculus'))
spc_carnivora = set(('canis_familiaris', 'felis_catus'))
spc_ev_ungl = set(('bos_taurus', 'sus_scrofa', 'ovis_aries')) # even toed ungulates
spc_rabbit = set(('oryctolagus_cuniculus',))
spc_horse = set(('equus_caballus',))


spc_map = {spc.split('_')[0].capitalize() + ' ' + spc.split('_')[1] : spc for spc in spc_all}
print(spc_map)


account = HostAccount('biodb01', 'benjamin_vernot', 'benjamin_vernot')
compara = Compara(spc_all, release=75, account=account)

print('SHOULD HAVE TWO HUMAN')
## shoudl have only two human regions
syntenic_regions = compara.get_syntenic_regions(species="Human", coord_name='1', start=25631543, end=25631594, align_method='EPO', align_clade='mammals')
for region in syntenic_regions:
    print(region)
    pass


print('SHOULD HAVE FOUR CHIMP')
## shoudl have four chimp regions
# 1:755186-755237
syntenic_regions = [r for r in compara.get_syntenic_regions(species="Human", coord_name='1', start=755186, end=755237, align_method='EPO', align_clade='mammals')]
print(len(syntenic_regions))
for region in syntenic_regions:
    print(region)
    pass


print('bad liftover?')
# 1:1080875-1080926
syntenic_regions = compara.get_syntenic_regions(species="Human", coord_name='1', start=1080875, end=1080926, align_method='EPO', align_clade='mammals')
print(list(syntenic_regions))



## window is +/- left and right (so default is 10 bases around target base)
def get_dists(aln_members, hum_seq, spc_set, is_anc = False, window = 5, hum_var_idx = None):
    
    current_species_set = set(spc_map[member.genome.species] for member in aln_members)
    min_indel_count = sys.maxsize
    min_snv_count = sys.maxsize
    if len(current_species_set.intersection(spc_set)) == 0:
        min_indel_count = 0
        min_snv_count = 0
        pass
    
    for member in aln_members:
        if spc_map[member.genome.species] not in spc_set: continue
        indel_count = sum(1 for i in range(len(hum_seq)) if member.aligned_seq[i] != hum_seq[i] and (hum_seq[i] == '-' or member.aligned_seq[i] == '-'))
        snv_count   = sum(1 for i in range(len(hum_seq)) if member.aligned_seq[i] != hum_seq[i] and hum_seq[i] != '-' and member.aligned_seq[i] != '-')
        print( 'INDEL', member.genome.species, indel_count, member.aligned_seq)
        min_indel_count = min(indel_count, min_indel_count)
        min_snv_count = min(snv_count, min_snv_count)
        pass
    
    return min_indel_count, min_snv_count


def report_target_base_state(spc, base, a1, a2, strand):
    primate_spc = ['Gorilla gorilla', 'Pan troglodytes', 'Pongo abelii', 'Homo sapiens']
    if spc in primate_spc and base != a1 and base != a2: return 'FAIL-TARGET %s : %s %s | S=%d' % (base, a1, a2, strand)
    if spc in primate_spc: return 'x | S=%d' % strand
    return ''

first_report = True
first_qc_report = True
idx = 0
for line in sys.stdin:
    idx += 1
    print(line.rstrip(), idx)

    chrom,pos = line.rstrip().split()[:2]
    pos = int(pos)

    if (chrom,pos) not in allelic_state:
        print( 'SITE NOT IN ALLELES FILE', chrom, pos)
        pass
    a1, a2 = allelic_state[(chrom,pos)][0], allelic_state[(chrom,pos)][1]


    if (chrom,pos) not in map_states:
        print( 'SITE NOT IN MAP/RPT FILE', chrom, pos)
        pass
    repeat_mask_bases, heng99_bases, heng99_target_overlap = map_states[(chrom,pos)]
    
    if (chrom,pos) not in ape_states:
        print( 'SITE NOT IN APES FILE', chrom, pos)
        pass
    hg19, pantro, panpan, gorgor, ponabe, rhemac, vind, altai, den = ape_states[(chrom,pos)]

    syntenic_regions = compara.get_syntenic_regions(species="Human", coord_name=chrom,
                                                    start=pos-offset, end=pos-offset+length,
                                                    align_method='EPO', align_clade='mammals')
    
    ## turn into a list
    aligned_pairs = [r for r in syntenic_regions]
    if len(aligned_pairs) == 0:
        print( 'SKIPPING REGION WITH NO HOMOLOGY %s:%d-%d' % (chrom,pos-offset,pos-offset+length))
        continue

    if len(aligned_pairs) > 1:
        print( 'SKIPPING REGION WITH %d HOMOLOGY BLOCKS %s:%d-%d' % (len(aligned_pairs),chrom,pos-offset,pos-offset+length))
        continue

    alignment = aligned_pairs[0]

    try:
        print( 'iter', idx, len(aligned_pairs), alignment.num_members, [a.region for a in alignment.members]) #alignment.getSpeciesSet())
    except AttributeError as err:
        print('SKIPPING DUE TO TRY/EXCEPT ERROR - not sure why this happens, but it\'s rare')
        print( err)
        continue
    
    #continue
    aln_members = [a for a in alignment.members if a.aligned_seq != None]
    all_species = [a.genome.species for a in aln_members]

    if len(set(all_species)) != len(all_species):
        print( 'EXTRA SPECIES')
        print( all_species)
        pass

    
    current_species = [spc_map[member.genome.species] for member in aln_members]
    current_species_set = set(current_species)
    hum_count = current_species.count('homo_sapiens')
    hum_idx = current_species.index('homo_sapiens')
    hum_seq = str(aln_members[hum_idx].aligned_seq)
    hum_seq_forward = (hum_seq if aln_members[hum_idx].location.strand == 1 else str(aln_members[hum_idx].aligned_seq.rc())).replace('-', '')

    #########
    ## find the base in the alignment that corresponds to our target base
    ## THIS DOESN'T REALLY WORK - ONLY WORKS IF THE STRAND IS 1
    ## IF STRAND OF HOMO IS -1, THEN SHOULD COUNT BACKWARDS IN THE ALIGNMENT!
    ## CURRENTLY ONLY USED FOR QC, AND IF OFFSET=1 [given as 0 on the command line, b/c I'm dumb and increment it immediately] AND LENGTH=1, SO WE ONLY FETCH THE TARGET BASE
    i = 0
    for aln_offset in range(len(hum_seq)):
        if hum_seq[aln_offset] != '-': i += 1
        if i == offset: break
        pass
    #########

    print( 'HUMSEQ', hum_seq.replace('-', ''), len(hum_seq.replace('-', '')))
    if len(hum_seq_forward) < length:
        print( 'SKIPPING REGION WITH SHORT HOMOLOGY BLOCK (%d bp) %s:%d-%d' % (len(hum_seq_forward),chrom,pos-offset,pos-offset+length))
        continue
        

    if len(current_species_set) < 2:
        print('Only one species, skipping')
        continue

    if hum_count > 1:
        print('TOO MANY HUMANS')
        continue

    # for a in [a for a in aln_members if a.aligned_seq != None]:
    x_hum_allele = None
    x_chimp_allele = None
    x_hum_strand = None
    x_primate_alleles = []
    x_nonprimate_alleles = []
    for a in aln_members:
        # print(dir(a))

        # if a.genome.Species != 'Sus scrofa': continue
        if a.aligned_seq == None:
            print( "skipping", a.genome.species)
            continue

        aln_seq_str = str(a.aligned_seq)
        
        if a.genome.species == 'Homo sapiens':
            print( aln_seq_str, a.genome.species, report_target_base_state(a.genome.species, aln_seq_str[aln_offset],
                                                                           a1, a2, a.location.strand))
            x_hum_allele = aln_seq_str[aln_offset]
            x_hum_strand = a.location.strand
        else:
            print( ''.join(aln_seq_str[i] if i == aln_offset or aln_seq_str[i] != hum_seq[i] else '.' for i in range(len(hum_seq))),
                   a.genome.species, report_target_base_state(a.genome.species, aln_seq_str[aln_offset],
                                                              a1, a2, a.location.strand))
            # print( ''.join(aln_seq_str[i] if i == aln_offset else ' ' for i in range(len(hum_seq))))
            if a.genome.species in ('Gorilla gorilla', 'Pan troglodytes', 'Pongo abelii'):
                x_primate_alleles += [aln_seq_str[aln_offset]]
            else:
                x_nonprimate_alleles += [aln_seq_str[aln_offset]]
                pass
            if a.genome.species == 'Pan troglodytes':
                x_chimp_allele = aln_seq_str[aln_offset]
                pass
            pass
        
        #print( a._cached['genome_db_id'], 'region' in a._cached)
        #print( a.region.genome, 'region' in a._cached, a._make_map_func) 
        pass

    if offset == 1 and length == 1:
        if first_qc_report == True:
            first_qc_report = False
            print('QC-ALLELE-COUNTS', 'hum_allele', 'a1', 'a2', 'hum_strand',
                  'apes_match_hg19',
                  'n_apes',
                  'ape_alleles',
                  'non_apes_match_hg19',
                  'n_non_apes',
                  'non_ape_alleles',
                  'non_apes_match_chimp')
            pass
        print('QC-ALLELE-COUNTS', x_hum_allele, a1, a2, x_hum_strand,
              sum([x == x_hum_allele and x != '-' for x in x_primate_alleles]),
              len(x_primate_alleles),
              '_'.join(sorted(set(x_primate_alleles))),
              sum([x == x_hum_allele and x != '-' for x in x_nonprimate_alleles]),
              len(x_nonprimate_alleles),
              '_'.join(sorted(set(x_nonprimate_alleles))),
              sum([x == x_chimp_allele and x != '-' for x in x_nonprimate_alleles]))
        pass

    

    idl_chimp, snv_chimp = get_dists(aln_members, hum_seq, spc_apes)
    idl_apes, snv_apes = get_dists(aln_members, hum_seq, spc_apes)
    idl_monkey, snv_monkey = get_dists(aln_members, hum_seq, spc_monkey)
    idl_rodents, snv_rodents = get_dists(aln_members, hum_seq, spc_rodents)
    idl_rabbit, snv_rabbit = get_dists(aln_members, hum_seq, spc_rabbit)
    idl_glires, snv_glires = get_dists(aln_members, hum_seq, spc_rodents.union(spc_rabbit))
    idl_carnivora, snv_carnivora = get_dists(aln_members, hum_seq, spc_carnivora)
    idl_ev_ungl, snv_ev_ungl = get_dists(aln_members, hum_seq, spc_ev_ungl)
    idl_horse, snv_horse = get_dists(aln_members, hum_seq, spc_horse)
    idl_ungl, snv_ungl = get_dists(aln_members, hum_seq, spc_ev_ungl.union(spc_horse))

    n_chimp = len(current_species_set.intersection(spc_chimp))
    n_apes, n_monkey = len(current_species_set.intersection(spc_apes)), len(current_species_set.intersection(spc_monkey))
    n_rodents, n_rabbit = len(current_species_set.intersection(spc_rodents)), len(current_species_set.intersection(spc_rabbit))
    n_glires = len(current_species_set.intersection(spc_rodents.union(spc_rabbit)))
    n_carnivora, n_ev_ungl = len(current_species_set.intersection(spc_carnivora)), len(current_species_set.intersection(spc_ev_ungl))
    n_horse, n_ungl = len(current_species_set.intersection(spc_horse)), len(current_species_set.intersection(spc_ev_ungl.union(spc_horse)))

    ## moved getting allelic state, etc, from here

    print( 'SEQDEBUG', hum_seq_forward)
    hum_seq_forward_a1 = hum_seq_forward[:offset-1] + a1 + hum_seq_forward[offset:]
    print( 'SEQDEBUG', hum_seq_forward_a1)
    hum_seq_forward_a2 = hum_seq_forward[:offset-1] + a2 + hum_seq_forward[offset:]
    print( 'SEQDEBUG', hum_seq_forward_a2)
    print( 'SEQDEBUG', (offset-1)*' ' + 'X', a1, a2, hum_seq_forward[offset-1] in allelic_state[(chrom,pos)])

    if not hum_seq_forward[offset-1] in allelic_state[(chrom,pos)]:
        print('Major issue with offset (possibly having to do with reversed sequence?)')
        print('Allele at offset site should match one of the two possible alleles', hum_seq_forward[offset-1], allelic_state[(chrom,pos)])
        sys.exit(-1)
        pass


    ## get number of clades w/ alignments (from n.apes, n.monkeys, n.glires, n.carnivora, n.ungl)
    n_5 = sum((n_apes>0, n_monkey>0, n_glires>0, n_carnivora>0, n_ungl>0))
    n_4 = sum((n_monkey>0, n_glires>0, n_carnivora>0, n_ungl>0))
    n_3 = sum((n_glires>0, n_carnivora>0, n_ungl>0))
    print( 'NCLADES', n_apes, n_monkey, n_glires, n_carnivora, n_ungl, 'N', n_5, n_4, n_3)

    ## get minimum burden over each relevant clade (5/4/3 w/ & w/o apes/monkeys)
    burdens = (snv_apes+idl_apes if n_apes>0 else sys.maxsize,
               snv_monkey+idl_monkey if n_monkey>0 else sys.maxsize,
               snv_glires+idl_glires if n_glires>0 else sys.maxsize,
               snv_carnivora+idl_carnivora if n_carnivora>0 else sys.maxsize,
               snv_ungl+idl_ungl if n_ungl>0 else sys.maxsize)
    b_5 = min(burdens)
    b_4 = min(burdens[1:])
    b_3 = min(burdens[2:])

    if first_report == True:
        first_report = False
        print(  'REPORT', 'chrom', 'pos', 'ref', 'a1', 'a2', \
                'region', 'aln_len', \
                'N', 'n_chimp', 'n_apes', 'n_monkey', 'n_rodents', 'n_rabbit', 'n_glires', 'n_carnivora', 'n_ev_ungl', 'n_horse', 'n_ungl', \
                'INDELS', 'idl_chimp', 'idl_apes', 'idl_monkey', 'idl_rodents', 'idl_rabbit', 'idl_glires', 'idl_carnivora', 'idl_ev_ungl', 'idl_horse', 'idl_ungl', \
                'SNVS', 'snv_chimp', 'snv_apes', 'snv_monkey', 'snv_rodents', 'snv_rabbit', 'snv_glires', 'snv_carnivora', 'snv_ev_ungl', 'snv_horse', 'snv_ungl', \
                'NCLADES', 'n_5', 'n_4', 'n_3', \
                'BURDEN', 'b_5', 'b_4', 'b_3', \
                'MAPRPT', 'repeat_mask_bases', 'heng99_bases', 'heng99_target_overlap', \
                'APES', 'hg19', 'pantro', 'panpan', 'gorgor', 'ponabe', 'rhemac', 'vind', 'altai', 'den', \
                'hum_seq_forward_a1', 'hum_seq_forward_a2')
        pass
    
    print(  'REPORT', chrom, pos, hum_seq_forward[offset-1], a1, a2, \
            "%s:%d-%d" % (chrom,pos-offset,pos-offset+length), len(hum_seq), \
            'N', n_chimp, n_apes, n_monkey, n_rodents, n_rabbit, n_glires, n_carnivora, n_ev_ungl, n_horse, n_ungl, \
            'INDELS', idl_chimp, idl_apes, idl_monkey, idl_rodents, idl_rabbit, idl_glires, idl_carnivora, idl_ev_ungl, idl_horse, idl_ungl, \
            'SNVS', snv_chimp, snv_apes, snv_monkey, snv_rodents, snv_rabbit, snv_glires, snv_carnivora, snv_ev_ungl, snv_horse, snv_ungl, \
            'NCLADES', n_5, n_4, n_3, \
            'BURDEN', b_5 if not b_5 == sys.maxsize else '.', b_4 if not b_4 == sys.maxsize else '.', b_3 if not b_3 == sys.maxsize else '.', \
            'MAPRPT', repeat_mask_bases, heng99_bases, heng99_target_overlap, \
            'APES', hg19, pantro, panpan, gorgor, ponabe, rhemac, vind, altai, den, \
            hum_seq_forward_a1, hum_seq_forward_a2)
    
    
    # print( alignment)
    # print( dir(alignment))
    # # aligned_regions = [m for m in alignment.members if m.region is not None]
    # aligned_regions = [m for m in alignment.members]
    # print( aligned_regions)
    # source_region, target_region = aligned_regions[:2]
    # print( dir(source_region.location))
    # print( source_region.location.coord_name, source_region.location.start, source_region.location.end)
    # print( target_region.location.coord_name, target_region.location.start, target_region.location.end)
    
    #print source_region.dirs()
    print()
    print()
    
    
print('SUCCESS')


if backbone_name == 'pUC-KanR_v1_mut1':
    # Create a dictionary mapping plasmid positions to backbone positions
    backbone_dict = {}
    for r in range(len(plasmid_range)):
        for i in range(plasmid_range[r][0], plasmid_range[r][1] + 1):
            backbone_position = i - plasmid_range[r][0] + backbone_range[r][0]
            backbone_dict[i] = backbone_position

    # Initialize lists to store results
    known_bb_snp = []
    known_bb_index = []

    # Iterate over listPosition and check if the value is within the backbone range
    for i, pos in enumerate(listPosition):
        if pos in backbone_dict:
            bb_index = backbone_dict[pos]
            if bb_index in [317, 875]:
                known_bb_index.append(i)
                known_bb_snp.append(listSNP[i])

    # Remove known BB SNPs from listPosition_snp and listSNP_snp
    for pos in sorted(known_bb_index, reverse=True):
        listSNP.pop(pos)
        listPosition.pop(pos)

    print(known_bb_snp)
    print(known_bb_index)

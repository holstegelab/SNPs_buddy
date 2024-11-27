import os
pj = os.path.join

input_dir = pj('/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED')

def get_regions(lrange):
    """Converts a region describer (tuple format) to a list of regions (string format).
    E.g. [('A', 1, 3, 1), ('A', 2, 4, 2)] -> ['A3H', 'A04']
    """
    res = []

    for component, level, splitnr, ploidy in lrange:
        if ploidy == 1:
            ploidy = 'H'
        else:
            ploidy = ''
        if level == 0:
            region = f'{component}{ploidy}'
        else:
            region = f'{component}{splitnr:0{level}d}{ploidy}'
        res.append(region)

    return res

level2_range_diploid_only = [('A', 2,x,2) for x in range(0,99)] + \
               [('X', 1,x,2) for x in range(0,5)] + \
               [('Y', 1, x,2) for x in range(0,2)]

level2_regions_diploid = get_regions(level2_range_diploid_only)
parts = level2_regions_diploid
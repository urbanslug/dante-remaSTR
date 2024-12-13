# data = json.loads(f)
# (sample, version, postfilter_params, motifs) = data
# (motif_id, nomenclatures, loci, phasing) = motifs[0]
# (locus_id, sequence, locus_nomenclatures, allele1, allele2, stats, spanning, flanking, graph_data) = loci[0]  # phasing[i] has the same structure
# (prediction, confidence, indels, mismatches, reads) = allele1  # allele2 has the same structure
# (confidence, indels, mismatches) = stats
#
# djlint dante_remastr_standalone_templates/population_report.html --reformat --format-css --format-js --close-void-tags --indent 2 --preserve-blank-lines

import json

data = json.load(open("./dante_remastr_standalone_templates/instance.json"))
print(data["sample"])
print(data["postfilter_params"])
for motif in data["motifs"]:
    print(motif)
    break

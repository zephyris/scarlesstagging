import sys
from scarlesstagging import ScarlessTagging

# python3 singleprimer.py gene_id terminus plasmid_system tag drug
# for example:
# python3 singleprimer.py Tb927.7.6580 both 2a mng bsr

# as in TriTrypDB
gene_ids = [sys.argv[1]]

#@markdown Tagging terminus
# "n", "c" or "both"
terminus = sys.argv[2]

# defines primer binding site sequences
# "ppot-compatible" or "2a"
plasmid_system = sys.argv[3]

if plasmid_system == "2a":
    # "msc" or "mng"
    tag = sys.argv[4]
    # "pac" or "bsr"
    drug = sys.argv[5]
else:
    tag = None
    drug = None

scarlesstagging = ScarlessTagging()
primer_results = []
for gene_id in gene_ids:
    result = scarlesstagging.design_primers(gene_id, terminus, plasmid_system, tag, drug)
    primer_results.append(result)

    for err in ["uerr", "derr"]:
        if result[err] is None:
            result[err] = "None"
        else:
            result[err] = ", ".join(result[err])
    for key in result:
        print("\t".join([str(x) for x in [key, result[key]]]))

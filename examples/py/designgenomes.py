from scarlesstagging import ScarlessTagging
import requests, os
import concurrent.futures
from tqdm import tqdm

# base url
version = "68"
base_url = "https://tritrypdb.org/common/downloads/release-"+version+"/"

outdir = "tritrypdb-genomes-v" + str(version)
if not os.path.exists(outdir):
    os.mkdir(outdir)

species = [
    "TbruceiTREU927",
    "TbruceiLister427",
    "TbruceigambienseDAL972",
    "TcongolenseIL3000",
    "TvivaxY486"
    "TbruceiEATRO1125",
    "TbruceiLister427_2018",
    "TcongolenseIL3000_2019",
    "TcongolenseTc1_148"
]

for specie in species:
    print("-"*len(specie))
    print(specie)
    print("-"*len(specie))

    gene_ids = []
    url = base_url + specie + "/gff/data/TriTrypDB-"+version+"_"+specie+".gff"
    response = requests.get(url)
    lines = response.text.split("\n")
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        if fields[2] == "protein_coding_gene":
            parameters = dict(item.split("=") for item in fields[8].split(";"))
            gene_ids.append(parameters["ID"])

    terminus = "both"

    settings = [
        {
            "plasmid_system": "ppot-compatible",
        },
        """
        {
            "plasmid_system": "prext2a-scarless",
            "tag": "mng",
            "drug": "bsr",
        },
        {
            "plasmid_system": "prext2a-scarless",
            "tag": "msc",
            "drug": "pac",
        },
        {
            "plasmid_system": "prext2a-scarless",
            "tag": "mng",
            "drug": "pac",
        },
        {
            "plasmid_system": "prext2a-scarless",
            "tag": "msc",
            "drug": "bsr",
        },
        """
    ]

    scarlesstagging = ScarlessTagging()

    def process_setting(setting):
        primer_results = []
        plasmid_system = setting["plasmid_system"]
        if plasmid_system == "prext2a-scarless":
            tag = setting["tag"]
            drug = setting["drug"]
        else:
            tag = None
            drug = None

        if plasmid_system == "ppot-compatible":
            outfile = outdir + "/" + specie + "_" + plasmid_system + ".txt"
        else:
            outfile = outdir + "/" + specie + "_" + plasmid_system + "-" + tag + "-" + drug + ".txt"
        
        # check if already predicted
        partial = False
        if os.path.exists(outfile):
            # read data
            predicted = open(outfile, "r").read().split("\n")
            # parse lines to get predicted gene IDs
            if len(predicted) > 0:
                partial = True
            predicted = [x.split("\t")[0] for x in predicted[1:] if x != "" and x != "###completed###"]
        # open output, appending
        file = open(outfile, "a")
        if partial:
            # write header 
            file.write("\t".join(["Gene ID", "UF", "UR", "USG", "UERR", "DF", "DR", "DSG", "DERR"])+"\r\n")
        # for all gene ids
        for gene_id in gene_ids:
            # unless already predicted
            if gene_id in predicted:
                continue
            # do the prediction
            result = scarlesstagging.design_primers(gene_id, terminus, plasmid_system, tag, drug)
            primer_results.append(result)
            # print a nicely formatted output to show progress
            for err in ["uerr", "derr"]:
                if result[err] is None:
                    result[err] = "None"
                else:
                    result[err] = ", ".join(result[err])
            for key in result:
                print("\t".join([str(x) for x in [key, result[key]]]))
            print("")
            # write output line
            file.write("\t".join([str(x) for x in [result["id"], result["uf"], result["ur"], result["usg"], result["uerr"], result["df"], result["dr"], result["dsg"], result["derr"]]])+"\r\n")
        if partial:
            file.write("###completed###" + "\r\n")
        file.close()
    
    """
    for setting in settings:
        process_setting(setting)
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for setting in settings:
            future = executor.submit(process_setting, setting)
            futures.append(future)

        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            pass
    
from scarlesstagging import ScarlessTagging
import requests
import concurrent.futures
from tqdm import tqdm

# base url
version = "68"
base_url = "https://tritrypdb.org/common/downloads/release-"+version+"/"

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
        {
            "plasmid_system": "2a",
            "tag": "mng",
            "drug": "bsr",
        },
        {
            "plasmid_system": "2a",
            "tag": "msc",
            "drug": "pac",
        },
        {
            "plasmid_system": "2a",
            "tag": "mng",
            "drug": "pac",
        },
        {
            "plasmid_system": "2a",
            "tag": "msc",
            "drug": "bsr",
        },
    ]

    scarlesstagging = ScarlessTagging()

    def process_setting(setting):
        primer_results = []
        plasmid_system = setting["plasmid_system"]
        if plasmid_system == "2a":
            tag = setting["tag"]
            drug = setting["drug"]
        else:
            tag = None
            drug = None

        if plasmid_system == "ppot-compatible":
            file = open(specie+"_"+plasmid_system+".txt", "w")
        else:
            file = open(specie+"_"+plasmid_system+"-"+tag+"-"+drug+".txt", "w")
        file.write("\t".join(["Gene ID", "UF", "UR", "USG", "UERR", "DF", "DR", "DSG", "DERR"])+"\r\n")
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
            print("")
            file.write("\t".join([str(x) for x in [result["id"], result["uf"], result["ur"], result["usg"], result["uerr"], result["df"], result["dr"], result["dsg"], result["derr"]]])+"\r\n")
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

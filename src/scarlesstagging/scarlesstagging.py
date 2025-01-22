
from urllib.request import urlopen
import math, time

class ScarlessTagging:
    def __init__(self):
        self.get_sequence_cache = {}
        self.minimum_homology_length = 30
        pass
  
    def sequence_request(self, id: str, upAnchor: str, upSign: str, upOffs: int, downAnchor: str, downSign: str, downOffs: int) -> str:
        """
        Fetch a sequence from TriTrypDB, behaviour as for sequence retrieval tool at https://tritrypdb.org/tritrypdb/app/fasta-tool/gene

        :param id: Gene id
        :param upAnchor: Upstream anchor: `"Start"` transcript start, `"CodeStart"` CDS start, `"CodeEnd"` CDS end, or `"End"` transcript end
        :param upSign: Upstream offset direction: `"plus"` (downstream) or `"minus"` (upstream)
        :param upOffs: Upstream offset distance in bases.
        :param downAnchor: Downstream anchor: options as for `upAnchor`
        :param downSign: Downstream offset direction: options as for `upSign`
        :param downOffs: Downstream offest distances in bases.
        :return: Sequence at the requested locus
        """
        # construct url
        url = "https://veupathdb.org/cgi-bin/geneSrt?project_id=EuPathDB"
        url += "&type=genomic" + "&downloadType=plain"
        url += "&ids=" + id
        url += "&upstreamAnchor=" + upAnchor
        url += "&upstreamSign=" + upSign + "&upstreamOffset=" + str(upOffs)
        url += "&downstreamAnchor=" + downAnchor
        url += "&downstreamSign=" + downSign + "&downstreamOffset=" + str(downOffs)
        # retrieve from cache, if cached
        if url in self.get_sequence_cache:
            return self.get_sequence_cache[url]
        # grab from url
        max_retries = 1000
        retry = 0
        while retry < max_retries:
            try:
                response = urlopen(url)
                string = response.read().decode(response.info().get_param("charset") or "utf-8-sig")
                result = "".join(string.splitlines()[1:])
                break
            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except:
                print("Failed to retrieve sequence, retrying. Attempt", retry, "of", max_retries)
                time.sleep(1)
                retry += 1
        if retry >= max_retries:
            raise Exception("Failed to retrieve sequence after", max_retries, "attempts")
        # cache response and return
        self.get_sequence_cache[url] = result
        return result

    def get_sequence(self, gene_id: str, feature: str, feature_end: str, indices: list):
        """
        Fetch a sequence from TriTrypDB, behaviour as for sequence retrieval tool at https://tritrypdb.org/tritrypdb/app/fasta-tool/gene

        :param id: Gene id
        :param feature: Reference feature: `"cds"` for coding sequence, `"transcript"` for transcript.
        :param end: Reference feature end: `"start"` for upstream end or `"end"` for downstream end.
        :param start_index:" Start of sequence relative to reference `feature` and `end`.
        :param downOffs: End of sequence relative to reference `feature` and `end`.
        :return: Sequence at the requested locus
        """
        # dict for translating features into request anchors
        anchor_translation = {
            "cds": {"start": "CodeStart", "end": "CodeEnd"},
            "transcript": {"start": "Start", "end": "End"}
        }
        # determine signs and offsets
        signs = [None, None]
        offs = [None, None]
        for index, value in enumerate(indices):
            if value >= 0:
                signs[index] = "plus"
                offs[index] = value
            else:
                signs[index] = "minus"
                offs[index] = -value
        offs[1] -= 1
        if signs[1] == "minus":
            offs[1] += 2
        # get sequence, handling broken behaviour of requests spanning 0
        seq = self.sequence_request(
            gene_id,
            anchor_translation[feature][feature_end], signs[0], offs[0],
            anchor_translation[feature][feature_end], signs[1], offs[1]
        )
        # check response for problems
        warnings = []
        request_name = gene_id+" "+feature+":"+str(indices[0])+"::"+str(indices[1])
        if "N" in seq.upper():
            warnings.append("Ns in " + request_name)
            seq = seq.replace("N", "")
        if len(seq) != indices[1] - indices[0]:
            warnings.append("Incorrect sequence length for "+request_name)
        # return result
        if len(warnings) > 0:
            return seq, warnings
        else:
            return seq, None

    def reverse_complement(self, seq) -> str:
        """
        Reverse complements a DNA/RNA sequence, preserving case. U replaced with T.

        :param seq: DNA or RNA sequence (without IUPAC degenerate bases)
        :return: Reverse complemented sequence
        """
        # replace U with T
        replace = {"U": "T"}
        for key in replace:
            seq = seq.replace(key.lower(), replace[key].lower())
        # define sequence complement
        complement = {"C": "G", "G": "C", "A": "T", "T": "A"}
        # update with lower case complement
        lower_complement = {}
        for key in complement:
            lower_complement.update({key.lower(): complement[key].lower()})
        complement.update(lower_complement)
        # generate reverse complement sequence
        rcseq = ""
        for c in range(len(seq)):
            rcseq = complement[seq[c]] + rcseq
        return rcseq

    def translate(self, seq: str) -> str:
        """
        self.translates a DNA sequence to protein sequence using global variable `codons`.

        :param seq: DNA sequence.
        :return: Protein sequence, in whatever translation convention used by `codons`.
        """
        # dna codon to single letter amino acid code dict
        codons = {
            'TAA': '*', 'TAG': '*', 'TGA': '*',                                     # *
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',                         # A
            'TGT': 'C', 'TGC': 'C',                                                 # C
            'GAT': 'D', 'GAC': 'D',                                                 # D
            'GAA': 'E', 'GAG': 'E',                                                 # E
            'TTT': 'F', 'TTC': 'F',                                                 # F
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',                         # G
            'CAT': 'H', 'CAC': 'H',                                                 # H
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I',                                     # I
            'AAA': 'K', 'AAG': 'K',                                                 # K
            'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', # L
            'ATG': 'M',                                                             # M
            'AAT': 'N', 'AAC': 'N',                                                 # N
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',                         # P
            'CAA': 'Q', 'CAG': 'Q',                                                 # Q
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', # R
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S', # S
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',                         # T
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',                         # V
            'TGG': 'W',                                                             # W
            'TAT': 'Y', 'TAC': 'Y'                                                  # Y
        }

        pseq = ""
        for codon_index in range(math.floor(len(seq) / 3)):
            pseq += codons[seq[codon_index * 3:codon_index * 3 + 3]]
        return pseq

    def find_mutatable_pam(self, seq: str):
        """
        Find PAM sites in a coding sequence which can be removed with a synonymous mutation.

        :param seq: DNA sequence of CDS, must be in frame with a whole number of codons (ie. length is a multiple of 3).
        :return: List of lists, recording pam position (start of NGG), direction (forward for NGG, reverse for CCN), mutation position and base for mutation
        """
        # list of dinucleotides which define pam
        pam_bases = ["GG", "CC"]

        # check length
        warnings = []
        if len(seq) % 3 != 0:
            warnings.append("Sequence for PAM search not a multiple of 3 in length (unexpected behaviour may occur)")
        # result list of mutatable pams
        result = []
        # make query sequence uppercase
        seq = seq.upper()
        pseq = self.translate(seq)
        # for each position in the sequence, if in a pam dinucleotide
        for i in range(len(seq)):
            for pam_offs in range(2):
                if seq[i:i+2] in pam_bases:
                    # check all mutations at that position
                    for r in ["A", "T", "C", "G"]:
                        if r != seq[i]:
                            mutseq = seq[:i] + r + seq[i + 1:]
                            if self.translate(mutseq) == pseq:
                            # determine pam position and direction from pam_offs and if GG vs. CC
                                if seq[i] == "G":
                                    pam_position = i - pam_offs - 1
                                    pam_cutpos = pam_position - 4
                                    pam_dir = "forward"
                                elif seq[i] == "C":
                                    pam_position = i - pam_offs + 3
                                    pam_cutpos = pam_position + 4
                                    pam_dir = "reverse"
                                result.append([pam_position, pam_dir, pam_cutpos, i, seq[i], r])
        if len(warnings) == 0:
            warnings = None
        return result, warnings

    # gene metadata check
    def check_gene_metadata(self, gene_id):
        url = "https://tritrypdb.org/tritrypdb/service/record-types/transcript/searches/GenesByText/reports/standard?text_search_organism=%5B%22Bodo%20saltans%20strain%20Lake%20Konstanz%22%2C%22Angomonas%20deanei%20strain%20Cavalho%20ATCC%20PRA-265%22%2C%22Blechomonas%20ayalai%20B08-376%22%2C%22Crithidia%20fasciculata%20strain%20Cf-Cl%22%2C%22Endotrypanum%20monterogeii%20strain%20LV88%22%2C%22Leishmania%20aethiopica%20L147%22%2C%22Leishmania%20amazonensis%20MHOM%2FBR%2F71973%2FM2269%22%2C%22Leishmania%20arabica%20strain%20LEM1108%22%2C%22Leishmania%20braziliensis%20MHOM%2FBR%2F75%2FM2903%22%2C%22Leishmania%20braziliensis%20MHOM%2FBR%2F75%2FM2904%22%2C%22Leishmania%20braziliensis%20MHOM%2FBR%2F75%2FM2904%202019%22%2C%22Leishmania%20donovani%20BPK282A1%22%2C%22Leishmania%20donovani%20CL-SL%22%2C%22Leishmania%20donovani%20HU3%22%2C%22Leishmania%20donovani%20strain%20LV9%22%2C%22Leishmania%20enriettii%20MCAV%2FBR%2F2001%2FCUR178%22%2C%22Leishmania%20enriettii%20strain%20LEM3045%22%2C%22Leishmania%20gerbilli%20strain%20LEM452%22%2C%22Leishmania%20infantum%20JPCM5%22%2C%22Leishmania%20major%20Friedlin%202021%22%2C%22Leishmania%20major%20strain%20Friedlin%22%2C%22Leishmania%20major%20strain%20LV39c5%22%2C%22Leishmania%20major%20strain%20SD%2075.1%22%2C%22Leishmania%20martiniquensis%20LEM2494%22%2C%22Leishmania%20martiniquensis%20MHOM%2FTH%2F2012%2FLSCM1%22%2C%22Leishmania%20mexicana%20MHOM%2FGT%2F2001%2FU1103%22%2C%22Leishmania%20orientalis%20MHOM%2FTH%2F2014%2FLSCM4%22%2C%22Leishmania%20panamensis%20MHOM%2FCOL%2F81%2FL13%22%2C%22Leishmania%20panamensis%20strain%20MHOM%2FPA%2F94%2FPSC-1%22%2C%22Leishmania%20sp.%20Ghana%20MHOM%2FGH%2F2012%2FGH5%22%2C%22Leishmania%20sp.%20Namibia%20MPRO%2FNA%2F1975%2F252%2FLV425%22%2C%22Leishmania%20tarentolae%20Parrot%20Tar%20II%202019%22%2C%22Leishmania%20tarentolae%20Parrot-TarII%22%2C%22Leishmania%20tropica%20L590%22%2C%22Leishmania%20turanica%20strain%20LEM423%22%2C%22Leptomonas%20pyrrhocoris%20H10%22%2C%22Leptomonas%20seymouri%20ATCC%2030220%22%2C%22Paratrypanosoma%20confusum%20CUL13%22%2C%22Porcisia%20hertigi%20MCOE%2FPA%2F1965%2FC119%22%2C%22Trypanosoma%20brucei%20EATRO1125%22%2C%22Trypanosoma%20brucei%20Lister%20strain%20427%22%2C%22Trypanosoma%20brucei%20Lister%20strain%20427%202018%22%2C%22Trypanosoma%20brucei%20brucei%20TREU927%22%2C%22Trypanosoma%20brucei%20gambiense%20DAL972%22%2C%22Trypanosoma%20congolense%20IL3000%22%2C%22Trypanosoma%20congolense%20IL3000%202019%22%2C%22Trypanosoma%20cruzi%20Berenice%22%2C%22Trypanosoma%20cruzi%20Brazil%20A4%22%2C%22Trypanosoma%20cruzi%20CL%20Brener%20Esmeraldo-like%22%2C%22Trypanosoma%20cruzi%20CL%20Brener%20Non-Esmeraldo-like%22%2C%22Trypanosoma%20cruzi%20Dm28c%202014%22%2C%22Trypanosoma%20cruzi%20Dm28c%202017%22%2C%22Trypanosoma%20cruzi%20Dm28c%202018%22%2C%22Trypanosoma%20cruzi%20Sylvio%20X10%2F1%22%2C%22Trypanosoma%20cruzi%20Sylvio%20X10%2F1-2012%22%2C%22Trypanosoma%20cruzi%20TCC%22%2C%22Trypanosoma%20cruzi%20Y%20C6%22%2C%22Trypanosoma%20cruzi%20marinkellei%20strain%20B7%22%2C%22Trypanosoma%20cruzi%20strain%20CL%22%2C%22Trypanosoma%20cruzi%20strain%20CL%20Brener%22%2C%22Trypanosoma%20cruzi%20strain%20G%22%2C%22Trypanosoma%20equiperdum%20OVI%22%2C%22Trypanosoma%20evansi%20strain%20STIB%20805%22%2C%22Trypanosoma%20grayi%20ANR4%22%2C%22Trypanosoma%20melophagium%20St.%20Kilda%22%2C%22Trypanosoma%20rangeli%20SC58%22%2C%22Trypanosoma%20theileri%20isolate%20Edinburgh%22%2C%22Trypanosoma%20vivax%20Y486%22%2C%22Eubodonida%22%2C%22Leishmania%20braziliensis%22%2C%22Leishmania%20donovani%22%2C%22Leishmania%20enriettii%22%2C%22Leishmania%20major%22%2C%22Leishmania%20martiniquensis%22%2C%22Leishmania%20panamensis%22%2C%22Leishmania%20tarentolae%22%2C%22Leptomonas%22%2C%22Trypanosoma%20brucei%22%2C%22Trypanosoma%20congolense%22%2C%22Trypanosoma%20cruzi%22%5D&text_expression="+gene_id+"&document_type=gene&text_fields=%5B%22primary_key%22%5D&reportConfig=%7B%22attributes%22%3A%5B%22primary_key%22%2C%22cds_length%22%2C%22is_pseudo%22%2C%22signalp_peptide%22%2C%22gene_type%22%5D%2C%22tables%22%3A%5B%5D%2C%22attributeFormat%22%3A%22text%22%7D"
        response = urlopen(url)
        string = response.read().decode(response.info().get_param("charset") or "utf-8-sig")
        import json
        data = json.loads(string)
        warnings = []
        warnings_n = []
        errors = ["Gene ID not found"]
        for record in data["records"]:
            if record["displayName"] == gene_id:
                errors = []
            if record["attributes"]["gene_type"] != "protein coding gene":
                warnings.append("Not a protein coding gene")
            if record["attributes"]["cds_length"] is None:
                warnings.append("No CDS length listed")
                record["attributes"]["cds_length"] = 0
            if int(record["attributes"]["cds_length"]) < self.minimum_homology_length * 2:
                warnings.append("Very short CDS may lead to misplaced sgRNAs")
            if record["attributes"]["signalp_peptide"] is not None:
                warnings_n.append("SignalP-predicted signal peptide "+str(record["attributes"]["signalp_peptide"])+" may prevent N terminal tagging")
            if record["attributes"]["is_pseudo"] != "No":
                warnings.append("Listed as a pseudogene")
        if len(warnings) == 0:
            warnings = None
        if len(warnings_n) == 0:
            warnings_n = None
        if len(errors) == 0:
            errors = None
        return errors, warnings, warnings_n
    
    def design_primers(self, gene_id: str, terminus: str, plasmid_system: str, tag: str, drug: str):
        tag = tag.lower()
        drug = drug.lower()

        # guide rna length (excluding pam ngg), bases
        guide_length = 20

        # minimum homology arm length, bases
        minimum_homology_length = 30

        # maximum primer length, bases
        maximum_primer_length = 100

        # primer binding sites
        binding_c_forward = None
        binding_c_reverse = None
        binding_n_forward = None
        binding_n_reverse = None
        if plasmid_system == "ppot-compatible":
            binding_n_forward = "gtataatgcagacctgctgc".lower()
            binding_n_reverse = "ctggatcaggatcgggtagt".lower()
            binding_c_forward = "ggttctggtagtggttccgg".lower()
            binding_c_reverse = "gcacaggtctctcaaattgg".lower()
            print("Warning: This primer design does not give scarless integration")
        elif plasmid_system == "prext2a-scarless":
            if tag == "mng":
                binding_n_reverse = "gtgcaagcgagcttggcg".lower()
                binding_c_forward = "ggtaccggttctggtagt".lower()
            if tag == "msc":
                binding_n_reverse = "gtgcaagcgagcttggcg".lower()
                binding_c_forward = "ggtaccggttctggtagt".lower()
            if drug == "bsr":
                binding_n_forward = "atgcctttgtctcaagaagaa".lower()
                binding_c_reverse = "ACCAATCAAGATCCCTTGGATTAA".lower()
            if drug == "pac":
                binding_n_forward = "ATGACTGAATACAAGCCAACG".lower()
                binding_c_reverse = "ACCAATCAAGATCCCTTGGATTAA".lower()
        else:
            raise ValueError("Plasmid system \" + plasmid_system + \" not recognised")
        if binding_n_forward is None or binding_n_reverse is None or binding_c_forward is None or binding_c_reverse is None:
            raise ValueError("Plasmid system \" + plasmid_system + \" not recognised")

        sgrna_t7 = "gaaattaatacgactcactatagg".lower()
        sgrna_hairpin = "gttttagagctagaaatagc".lower()

        verbose = False

        #@title Find tagging primers

        #@markdown Exhaustively scans the start (for N terminal tagging) or end (for C terminal tagging) of the open reading frame for PAM sites which can be mutated to remove the PAM without changing the coding sequence.

        #@markdown Up to six types of primers are designed:

        def print_cds_pam_result(seq, pam_results):
            print("PAM search area ", seq)
            print("Protein sequence", self.translate(seq))
            print("Results:", "       ", "("+", ".join(["PAM posn.", "PAM dir.", "Cut posn.", "Mut. posn.", "Orig. base", "Mut. base"])+")")
            for result in pam_results:
                print("", ", ".join([str(x) for x in result]))

        primer_result = {"id": gene_id, "uf": None, "ur": None, "usg": None, "uerr": None, "df": None, "dr": None, "dsg": None, "derr": None}
        gene_errors, gene_warnings, gene_warnings_n = self.check_gene_metadata(gene_id)
        if gene_errors is not None:
            primer_result["uerr"] = gene_errors
            primer_result["derr"] = gene_errors
        else:
            if terminus == "n" or terminus=="both":
                #@markdown For N terminal tagging, gives an upstream forward (UF), upstream reverse (UR) and upstream sgRNA (USG) primer
                # find pams in cds, start of protein coding sequence
                if gene_warnings is not None and gene_warnings_n is not None:
                    warnings = gene_warnings.copy() + gene_warnings_n.copy()
                else:
                    warnings = []
                pam_search_range = 3 * math.floor((maximum_primer_length - len(binding_n_reverse) - minimum_homology_length) / 3)
                seq, warning = self.get_sequence(gene_id, "cds", "start", [3, 3 + pam_search_range])
                if warning: warnings += warning
                pam_results, warning = self.find_mutatable_pam(seq)
                if warning: warnings += warning
                if verbose: print_cds_pam_result(seq, pam_results)
                if len(pam_results) == 0:
                    # if no sites found, fail!
                    primer_result["uerr"] = ["No mutatable PAM sites found near start of CDS"]
                else:
                    # first in list, closest to n terminus
                    pam_result = pam_results[0]
                    # upstream forward primer
                    upstream_forward_sequence, warning = self.get_sequence(gene_id, "cds", "start", [-minimum_homology_length, 0]) # Removed offset of -1
                    if warning: warnings += warning
                    upstream_forward = upstream_forward_sequence + binding_n_forward
                    if verbose: print("Upstream forward:    ", upstream_forward)
                    primer_result["uf"] = upstream_forward
                    # upstream reverse primer
                    reverse_homology, warning = self.get_sequence(gene_id, "cds", "start", [3, 3 + minimum_homology_length + pam_result[2]])
                    if warning: warnings += warning
                    reverse_homology = reverse_homology[:pam_result[3]] + pam_result[5].lower() + reverse_homology[pam_result[3] + 1:]
                    upstream_reverse = self.reverse_complement(binding_n_reverse + reverse_homology)
                    if verbose: print("Upstream reverse:    ", upstream_reverse)
                    primer_result["ur"] = upstream_reverse
                    # upstream sgrna primer
                    if pam_result[1] == "forward":
                        upstream_sgrna_sequence, warning = self.get_sequence(gene_id, "cds", "start", [3 + pam_result[0] - guide_length, 3 + pam_result[0]])
                        if warning: warnings += warning
                        upstream_sgrna = sgrna_t7 + upstream_sgrna_sequence + sgrna_hairpin
                    elif pam_result[1] == "reverse":
                        upstream_sgrna_sequence, warning = self.get_sequence(gene_id, "cds", "start", [3 + pam_result[0], 3 + pam_result[0] + guide_length])
                        if warning: warnings += warning
                        upstream_sgrna = sgrna_t7 + self.reverse_complement(upstream_sgrna_sequence) + sgrna_hairpin
                    if verbose: print("Upstream sgRNA:        ", upstream_sgrna)
                    primer_result["usg"] = upstream_sgrna
                    if len(warnings) > 0:
                        primer_result["uerr"] = warnings
            if terminus == "c" or terminus=="both":
                #@markdown For C terminal tagging, gives a downstream forward (DF), downstream reverse (DR) and downstream sgRNA (DSG) primer
                # find pam in cds, end of protein coding sequence, minus stop codon
                if gene_warnings is not None:
                    warnings = gene_warnings.copy()
                else:
                    warnings = []
                pam_search_range = 3 * math.floor((maximum_primer_length - len(binding_c_forward) - minimum_homology_length) / 3)
                seq, warning = self.get_sequence(gene_id, "cds", "end", [-pam_search_range - 3, -3])
                if warning: warnings += warning
                pam_results, warning = self.find_mutatable_pam(seq)
                if warning: warnings += warning
                if verbose: print_cds_pam_result(seq, pam_results)
                if len(pam_results) == 0:
                    #if no sites found, fail!
                    primer_result["derr"] = ["No mutatable PAM sites found near end of CDS"]
                else:
                    # last in list, closest to c terminus
                    pam_result = pam_results[-1]
                    # pam sites were found in forward direction, so adjust with pam_search_range
                    # downstream forward primer
                    forward_homology, warning = self.get_sequence(gene_id, "cds", "end", [-minimum_homology_length - (pam_search_range - pam_result[2]) - 3 + 1, - 3 + 1]) # Added offset of +1 to correct for gene end indexing
                    if warning: warnings += warning
                    forward_homology = forward_homology[:-pam_search_range + pam_result[3]] + pam_result[5].lower() + forward_homology[-pam_search_range + pam_result[3] + 1:] # TODO? Is an offset of +1 also required here? Is PAM site correct in position?
                    downstream_forward = forward_homology + binding_c_forward
                    if verbose: print("Downstream forward:", downstream_forward)
                    primer_result["df"] = downstream_forward
                    # downstream reverse primer
                    downstream_reverse_sequence, warning = self.get_sequence(gene_id, "cds", "end", [0 + 1, minimum_homology_length + 1])
                    if warning: warnings += warning
                    downstream_reverse = self.reverse_complement(binding_c_reverse + downstream_reverse_sequence)
                    if verbose: print("Downstream reverse:", downstream_reverse)
                    primer_result["dr"] = downstream_reverse
                    # downstream sgrna primer
                    if pam_result[1] == "forward":
                        downstream_sgrna_sequence, warning = self.get_sequence(gene_id, "cds", "end", [-pam_search_range + pam_result[0] - guide_length - 3 + 1, -pam_search_range + pam_result[0] - 3 + 1])
                        if warning: warnings += warning
                        downstream_sgrna = sgrna_t7 + downstream_sgrna_sequence + sgrna_hairpin
                    elif pam_result[1] == "reverse":
                        downstream_sgrna_sequence, warning = self.get_sequence(gene_id, "cds", "end", [-pam_search_range + pam_result[0] - 3 + 1, -pam_search_range + pam_result[0] + guide_length + 1]) # Added offset of +1 to correct for gene end indexing
                        if warning: warnings += warning
                        downstream_sgrna = sgrna_t7 + self.reverse_complement(downstream_sgrna_sequence) + sgrna_hairpin
                    if verbose: print("Downstream sgRNA:    ", downstream_sgrna)
                    primer_result["dsg"] = downstream_sgrna
                    if len(warnings) > 0:
                        primer_result["derr"] = warnings
        return primer_result

# Scarless Tagging
This is a primer design tool for endogenous tagging in trypanosomatid parasites retaining unmodified endogenous mRNA untranslated regions (UTRs).
Primers are designed for the pRExT2A plasmid system, using either pPOT-compatible (which leave a small scar) or drug and tag-specific primers for scarless tagging.

To design primers for your genes of interest, go to [this collaboratory notebook](https://colab.research.google.com/github/zephyris/scarlesstagging/blob/main/examples/ipynb/ScarlessTagging.ipynb) and follow the instructions. This runs in your browser - so no installation to worry about!

The primer design tool is written as a Python module which can also be installed for more advanced usage. Install using `python3 -m pip install git+https://github.com/zephyris/scarlesstagging`. You can find some example code in `examples/py/`.

This primer design tool fetches sequence data from [TriTrypDB.org](https://tritrypdb.org/) and depends on their continued funding and operation to work.

The full description of how this tagging system works, how the primer design strategy works and validation of the system is available on BioRxiv: [Carbajo et al. 2025](https://doi.org/10.1098/rsob.240334).
If you use this tool for your work then please cite this paper.

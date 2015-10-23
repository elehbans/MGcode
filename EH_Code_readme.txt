The following code adds on to the work by the Salis lab regarding the analysis and design of Ribosome Binding Sites:

RBS_CalculatorEH.py                   * a modified version to report more variables and data to other programs
EH_parser.py                          * Code to parse a FASTA formate genome using a corresponding genbank file for all genes (-35 to +35 of start codon)
Matrix_Run_fullgenome.py              * Takes the resulting file from EH_parser.py (a txt in FASTA format) and runs all genes through the RBS_CalculatorEH.py, outputs csv
Matrix_design_MTM.py                  * Allows user to design new (-35 to +35) sequence using command-line arguments and outputs tracker file
Matrix_batch_design.py                * Allows user to do batches of designs (e.g. 10 versions) for a number of input -35/+35 sections
Matrix_design_module.py               * A converted version of Matrix_design_MTM.py that is used by Matrix_batch_design.py
spacing_ids.csv                       * A necessary CSV for reporting the spacing of RBS in nt instead of dG_spacing
Chlamy_codons.csv                     * A csv with codon frequencies for Chlamydomonas Reinhardtii

# saier_lab

## [Python Program](https://github.com/CobyLsy/saier_lab/blob/main/mmseqs_hmmtop.py)
- Program that runs `mmseqs search`, `mmseqs convertalis`, and `hmmtop` on designated files; also pulls substrate data from tcdb.org.
- Result files from `convertalis`,`hmmtop`, and substrate data will be organized and stored into three dictionary objects (see **Sample Run** for example)
- Program asks for 3 input files (2 for both `search` and `convertalis` and 1 for `hmmtop`), and 3 names for output files (1 for both `search` and `convertalis`, 1 for `hmmtop`, and 1 for substrate data)

## [Sample Run](https://github.com/CobyLsy/saier_lab/blob/main/Sample-Run)
- 3 lines were the output when `print()` is called on 3 dictionaries, in which results of `convertalis`, `hmmtop`, and substrate data were stored
- Files and file names used in the run that produces these results: 1. QueryDB: `proteomeDB` 2. TargetDB: `tcdbDB` 3. Result Name: `pyDB` 4. Fasta file: `sample.faa` (from first 100 lines of `tcdbDB`) 5. hmmtop Result Name: `sample_hmm.hmmtop` 6. Substrate data file Name: `substrate.tsv`

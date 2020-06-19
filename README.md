# rpSelenzyme

* Docker Image: [brsynth/rpselenzyme-standalone](https://hub.docker.com/r/brsynth/rpselenzyme-standalone)

Tool that takes for input a tar.xz with a collection of rpSBML or single rpSBML, scans for the reaction rules and makes a REST request to [Selenzyme](Selenzyme.synbiochem.co.uk) and finds the enzymatic sequences (through Uniprot ID's) of all the reactions in heterologous pathways of rpSBML files.

## Input

Required:
* **-input**: (string) Path to the input file
* **-input_format**: (string) Valid options: sbml, tar. Format of the input file
* **-taxonomy_input**: (integer) Host taxomonomy ID of the organism the enzymes will be expressed in
* **-taxonomy_format**: (string) Valid options: json, integer. Format of the taxonomy input

Advanced Options:
* **-pathway_id**: (string, default: rp_pathway)
* **-num_results**: (integer, default: 10) Number of uniprot ID's to assign per reaction 
* **-direction**: (integer, default: 0) Keep only direction in the database (MetaNetX) or compare both forward and reverse
* **-noMSA**: (boolean, default: True) Do not perform multiple sequence alignement
* **-fp**: (string, default: RDK) Chemical fingerprint used in the similarity algorithm
* **-rxntype**: (string, default: smarts) Type of reaction rules that are passed to Selenzyme
* **-min_aa_length**: (integer, default: 100) Minimal amino acid length of the enzymatic sequences

## Ouput

* **-output**: (string) Path to the output file

## Dependencies

* Base Docker Image: [brsynth/rpBase](https://hub.docker.com/r/brsynth/rpbase)

## Installing 

WARNING: This tool requires a dataset that is not available to the public. To gain access to this data please contact the author of the tool.

To compile the docker run the following command:

```
docker build -t brsynth/rpselenzyme-standalone:dev .
```

## Running the test

To test the execution of the tool, untar the test.tar.xz and execute the following command:

```
python run.py -input test/test_rpThermo.tar -input_format tar -taxonomy_format string -taxonomy_input 83333 -min_aa_length 100 -output test/test_rpSelenzyme.tar
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

v0.1

## Authors

* **Pablo Carbonell**
* Melchior du Lac

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

### How to cite rpSelenzyme?

Carbonell, Pablo, et al. "Selenzyme: Enzyme selection tool for pathway design." Bioinformatics 34.12 (2018): 2153-2154.

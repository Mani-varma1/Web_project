# POPGen: A population genomics web application

## Link: 

## Description

The aim of this web application is to provide a proof-of-concept for an online application that can aid in population genomic analysis.

The data is based on The 1000 Genomes Project phase 3 v5b. RS IDs are based on the current build for dbSNP (build 155). 

Currently, this application allows the calculation of windowed statistics for nucleotide diversity, haplotype diversity, Tajima's D, and average Hudson Fst. Statistics are calculated based on user choice and use scitkit-allel for most of the functions. A link to their GitHub can be found below:

https://github.com/cggh/scikit-allel

Window size and step size are also user submitted. Overall statistics are provided per population and for the overall sample set (total populations selected). 

Plots are produced for each windowed statistic calculated using plotly. 

Data can be also downloaded for overall statistics and query results.

This application was developed as part of the QMUL MSc Bioinformatics Software Development Group Project. 


## Setting Up Locally
POPgen can run locally using: 

<code>git clone https://github.com/Mani-varma1/Web_project</code> 

<code>cd Flask</code>   

Virtual environment can be set up using:  

<code>py -3 -m venv .venv</code>

<code>.venv/Scripts/activate</code>  

To install necessary packages: 

<code>pip install -r packages.txt</code>  

A script for populating the database has been provided for those that wish to set up a local database and wish to explore the dataset. The files required for generating the database have also been provided in the form of a zip file with an empty database. All that needs to be done is to run the database_populating.py script in the folder it is located in. This will extract the files to a new data folder and automatically populate the database. The script takes roughly 30 minutes to run.  

<code>python database_populating.py</code>

Finally, to run the application: 

<code>python application.py</code>  

## Search Parameters

Searches can be made using genomic positions, RS ID, and gene names. Gene names have been limited to 3 for sake of computational requirements. 

Multiple RS IDs and gene names can be passed as comma-separated. Gene names will need to be given using HGNC conventions (case-sensitive). 


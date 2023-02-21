# LC-MS Target Compound Identification

This repository contains a script that identifies target compounds present in a dataset of LC-MS peaks and generates a database of identified compounds. The script allows you to analyze this database and visualize the results using an interactive Datasette webfront. 

## Requirements 

To use this script, you need to have the following software installed:

- Docker
- Python

## How to run

To get started, clone or download this repository. Next, install all Python dependencies using the command:
```
pip install -r requirements.txt
``` 

# Usage 

### Generating the database of identified compounds

To analyze a database of target compounds identified on the LC-MS data use the command:
```
python lcms_identification.py generate -l <peaklist.csv> -d <database.db> 
```
where ``<peaklist.csv>`` is a CSV file of experimental data containing (at least) the following columns for each peak in the LC-MS:
  - ``mz`` - Mass to charge ratio of the peak
  - ``rt`` - Retention time of the peak
  - ``intensity`` - Peak intensity

and ``<database.db>`` is a SQL database of target/suspect compounds containing (at least) a table ``compoundlist`` with (at least) the following columns:

  - ``mass_to_charge_ratio`` - Mass to charge ratio of the compound
  - ``retention_time`` - Retention time of the compound
  - ``retention_time_tolerance`` - Tolerance value for matching the compound's retention time

The script will identify the compounds present in ``<peaklist.csv>`` and generate a ``results_database.db`` database containing information on the identified compounds. Additionally, the script will open an interactive Datasette webfront to visualize the results.  

### Connecting to the Webfront

To directly connect to the Datasette webfront containing the results of the last dataset analyzed (without repeating the analysis), use the command:
```
python lcms_identification.py connect
```

### Downloading as CSV file

To download a CSV file containing the results of the results contained in the Datasette webfront, run the command:
```
python lcms_identification.py download
```

## Additional information 

To get more information about the script and its options, run the command:
```
python lcms_identification.py --help
```

## License 

This project is licensed under the MIT License.
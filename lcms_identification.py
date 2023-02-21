import numpy as np 
import sqlite3 as sql 
import pandas as pd 
import subprocess
import webbrowser
import argparse 
from time import sleep 

# Shortcut for readability
run = lambda command: subprocess.run(command.split())

# =======================================================================================================
def MSpeak_target_compound_identification(peaklist_csv_file, database_db_file, default_mass_tolerance=0.002, default_retime_tolerance=0.5):

    """
    Analyze MS peaks and identify target compounds from a database.

    Parameters
    ----------
    peaklist_csv_file :: str 
        Path to a CSV file containing the MS peak list.
    database_db_file :: str 
        Path to an SQL database file containing the target compound database.
    default_mass_tolerance
        Default mass-to-charge ratio tolerance (Da) for assigning the LC-MS peaks. Defaults to 0.002.
    default_retime_tolerance
        Default mass-to-charge ratio tolerance (Da) for assigning the LC-MS peaks. Defaults to 0.002.
    """

    def distance_map(quantity):
        return np.abs(np.atleast_2d(experimental[quantity]).T - np.atleast_2d(database[quantity]))

    def connection_map_based_on(quantity):
        return distance_map(quantity) < tolerance[quantity]

    def peaks_ids(target_id):
        return peak_target_matches[np.where(peak_target_matches[:,1]==target_id)[0]].flatten()[0]

    def target_intensities(target):
        return experimental['intensities'][peaks_ids(target)]

    def average_over_peaks(value, target):
        intensities = target_intensities(target)
        return np.sum(value*intensities/np.sum(intensities))

    def experimental_value_of_target(quantity,target_id):
        return average_over_peaks(experimental[quantity][peaks_ids(target_id)], target_id)

    def experimental_error_of_target(quantity,target_id):
        experimental[quantity][peaks_ids(target_id)]
        errors = np.abs(experimental[quantity][peaks_ids(target_id)] - database[quantity][target_id])
        if quantity=='m/z':
            errors = errors*1e6/database['m/z'][target_id]
        return average_over_peaks(errors, target_id)
      
    # Connect and retrieve target compound database
    con = sql.connect(database_db_file)
    targets_database = pd.read_sql_query("SELECT * FROM compoundlist",con)

    # Load the list of MS-peaks
    peaklist = pd.read_csv(peaklist_csv_file)

    # Data preparation
    # -----------------------------------------------------------------------------------------------------------

    database_labels = ['mass_to_charge_ratio','retention_time','retention_time_tolerance']
    newlabels = ['m/z','retime','retime_tolerance']
    database = {f'{newlabel}' : targets_database[label].to_numpy() for newlabel,label in zip(newlabels,database_labels) }

    # Convert None to np.nan
    database['retime_tolerance'] = database['retime_tolerance'].astype(float)

    # Extract ndarrays from the peaklist dataframe
    csv_labels = ['mz','rt','intensity']
    newlabels = ['m/z','retime','intensities']
    experimental = { f'{newlabel}' : peaklist[label].to_numpy() for newlabel,label in zip(newlabels,csv_labels) }

    # Define the tolerances for each quantity
    tolerance = {
        'm/z' : default_mass_tolerance,
        'retime' : database['retime_tolerance']
    }
    # Set the tolerance for retention time to the default value for any targets where it is NaN
    tolerance['retime'][np.isnan(database['retime_tolerance'])] = default_retime_tolerance


    # Analysis
    # -----------------------------------------------------------------------------------------------------------

    # Create a connection map for each quantity
    connection_map  = {label : connection_map_based_on(label) for label in ['m/z','retime'] }

    # Include those targets whose retention times are not known
    connection_map['retime'] = connection_map['retime'] | np.isnan(distance_map('retime'))
    
    # Find the peaks that match to the targets
    peak_target_matches = np.array( np.where( connection_map['m/z'] & connection_map['retime'] ) ).T

    # Find the targets that have been identified
    identified_targets, peaks_per_target = np.unique(peak_target_matches[:,1], return_counts=True)

    # Calculate the total intensity for each identified target and normalize it to a value in parts-per-million
    total_intensities = [np.sum(target_intensities(target_id)) for target_id in identified_targets]
    total_intensities = total_intensities/np.sum(experimental['intensities'])*1e6

    # Calculate the experimental values and errors for each identified compound
    compounds_values = { label: [experimental_value_of_target(label,target_id) for target_id in identified_targets] for label in ['m/z','retime'] }
    compounds_errors = { label: [experimental_error_of_target(label,target_id) for target_id in identified_targets] for label in ['m/z','retime'] }

    # Print the number of identified target compounds
    print(f'Identified {len(identified_targets)} target compounds.')

    # Output
    # -----------------------------------------------------------------------------------------------------------

    # Create a dataframe with identified targets data
    identified_targets_database = pd.DataFrame({
        'Compound ID' : targets_database.compound_id.loc[identified_targets],
        'Compound name' : targets_database.compound[identified_targets],
        'Total intensity (ppm)' : np.round(total_intensities).astype(int),
        'm/z (Da)' : np.round(compounds_values['m/z'],2),
        'm/z error (ppm)' :  np.round(compounds_errors['m/z']).astype(int),
        'Ret. time (min)' :  np.round(compounds_values['retime'],2),
        'Ret. time error (min)' :  np.round(compounds_errors['retime'],2),
        'Peaks within tolerance' :  peaks_per_target
    })
    # Sort identified targets dataframe by Total intensity (ppm) in descending order
    identified_targets_database = identified_targets_database.sort_values('Total intensity (ppm)',ascending=False)
    # Reset index of identified targets dataframe
    identified_targets_database = identified_targets_database.reset_index()

    # Define output name and connect to the database
    connection_results = sql.connect(f'results_database.db') 
    # Export identified compounds dataframe to the SQLite database
    identified_targets_database.to_sql('Compounds identified by LC-MS', connection_results, if_exists='replace')

    # Build Docker image
    run('docker build -f .\Dockerfile . -t datasette-web-front')    
# =======================================================================================================

# =======================================================================================================
def serve_datasette_webfront(port):
    """
    Starts a Datasette web front server in a Docker container on the specified port.

    Parameters
    ----------
    port :: int
        The port number to use for the web front server.
    """
    # Stop container service if still running from a previous execution
    run('docker stop datasette-web-front-service')
    # Replace existing container (if exists)
    run('docker rm datasette-web-front-service')
    # Create new container and start the web front server 
    run(f'docker run --name datasette-web-front-service -d -p {port}:8080 datasette-web-front')

    # Wait for some seconds for the webfront to be accessible
    sleep(3)
# =======================================================================================================

# =======================================================================================================
def connect_to_datasette_webfront(port):
    """
    Opens a web browser to a specific Datasette URL with predefined options.

    Parameters
    ----------
    port :: int
        The port number of the running Datasette web front server.
    """
    # Start the webfront
    serve_datasette_webfront(port) 

    # Set base URL for the database and define Datasette options
    base_url = f"http://localhost:{port}/results_database/Compounds+identified+by+LC-MS"
    default_sorting ='_sort_desc=Total+intensity+(ppm)'
    hide_index_columns = '&_nocol=level_0&_nocol=index'
    default_plot = '#g.mark=bar&g.x_column=Compound name&g.x_type=ordinal&g.y_column=Total intensity (ppm)&g.y_type=quantitative'

    # Create full URL with all the options defined and open it in a web browser
    url = base_url + '?' + default_sorting + hide_index_columns + default_plot
    webbrowser.open(url, new=0, autoraise=True)
# =======================================================================================================


# =======================================================================================================
def download_from_datasette_webfront(port):
    """
    Downloads a CSV file from a Datasette web front server to the local machine.

    Parameters
    ----------
    port :: int
        The port number of the running Datasette web front server.
    """
    # Start the webfront
    serve_datasette_webfront(port) 

    # Set URL for downloading the CSV from Datasette
    download_url = f"http://localhost:{port}/results_database/Compounds+identified+by+LC-MS.csv"

    # Download
    df = pd.read_csv(download_url)   
    df.to_csv('Compounds_identified_by_LC-MS.csv')
# =======================================================================================================


if __name__ == "__main__":

    # Create parser for input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("action", type=str, choices=['generate', 'connect', 'download'], 
                        help="Action to be performed:" + 
                                '\n\t generate - Generate a database of identified compounds and open a Datasette webfront server for visualization.'+
                                '\n\t open - Open the Datasette webfront server for visualization.'+
                                '\n\t download - Download a CSV file with the identified compounds database data from the Datasette Docker container.')
    parser.add_argument("-l", "--peaklist-file", type=str, help="Path to a CSV file containing the MS peak list. Required for action='generate'.")
    parser.add_argument("-d", "--database-file", type=str, help="Path to an SQL database file containing the target compound database. Required for action='generate'.")
    parser.add_argument("-p", "--port", type=int, default=8080, nargs='?', help="Port on which to open connection to Datasette webfront. Defaults to port 8080.")
    parser.add_argument("-m", "--mass-tolerance", type=float, default=0.002, nargs='?', help="Default mass-to-charge ratio tolerance (Da) for assigning the LC-MS peaks. Defaults to 0.002.")
    parser.add_argument("-t", "--time-tolerance", type=float, default=0.5, nargs='?', help="Default retention time tolerance (minutes) for assigning the LC-MS peaks. Defaults to 0.5min.")
    args = parser.parse_args()

    # For the 'generate' action, the filenames are required
    if args.action == 'generate' and not (args.peaklist_file and args.database_file):
        raise SyntaxError('To generate the webfront you must specify the --peaklist-file and --database-file arguments.')

    # Analyze the data and create the Datasette webfront image
    if args.action == 'generate':
        # Launch analysis
        MSpeak_target_compound_identification(args.peaklist_file,args.database_file, default_retime_tolerance=args.mass_tolerance, default_mass_tolerance=args.time_tolerance)

    # Connect to or download from the Datasette webfront
    if args.action == 'generate' or args.action == 'connect' :
        # Start web-front server with Datasette 
        connect_to_datasette_webfront(args.port)
    else: 
        # Download the data directly from datasette
        download_from_datasette_webfront(args.port)

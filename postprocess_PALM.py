from control_post_processor import ControlPostProcessor
import argparse
import os
import json

'''
-------------------------------------------------------------------
    Command-line interface for the PALM Ensight post-processor.    
-------------------------------------------------------------------
    Arguments
-------------------------------------------------------------------

folder, str
    Absolute path to the folder containing the fluent results.

parameters, str
    Json file containing the relevant parameters

options, str
    Json file containing the post-processing options

-------------------------------------------------------------------
    Usage
-------------------------------------------------------------------

Use options specified by a .json file
>>> uv run python -u postprocess_PALM.py "path/to/fluent/data" "parameters.json" -o "options.json"

Use options specified by user at runtime
>>> uv run python -u postprocess_PALM.py "path/to/fluent/data" "parameters.json"

See help menu
>>> uv run python postprocess_PALM.py -h

-------------------------------------------------------------------
'''


def print_postprocessing_setup(folder, options, parameters) -> None:
    '''Prints the command line options specified by the user prior to running.'''

    print('Executing Ensight post-processor with:')
    print(f'Data Folder = "{folder}"')

    print('Parameters of run: ')
    for item in parameters.items():
        print(f'{item[0]} = {item[1]}')

    if options:
        print('Post-processing options selected: ')
        for item in options.items():
            print(f'{item[0]} = {item[1]}')
    else:
        print('Options will be selected at runtime.')

    print('Starting post-processing:')


def main():
    '''Command line argument parser for PALM postprocessing'''

    parser = argparse.ArgumentParser(
        description = 'Post-process a PALM fluent run, from a predefined set of options.'
    )

    parser.add_argument('folder', help='Absolute path to the folder containing the fluent results.')
    parser.add_argument('parameters', help='json file containing the relevant parameters')
    parser.add_argument('-o', '--options', help='json file containing the post-processing options', metavar='*.json')

    args = parser.parse_args()
    

    if os.path.isdir(args.folder):
        folder_to_read = os.path.abspath(args.folder)
        print(f'The folder "{args.folder}" was specified.')

    else:
        raise FileNotFoundError(f'ERROR: The folder specified "{args.folder}" is not valid.')


    if os.path.isfile(args.parameters) and '.json' in args.parameters:

        parameters_file_name = args.parameters
        print(f'The parameter file "{args.parameters}" was specified.') 

        with open(parameters_file_name) as parameters_file:
            parameters = json.load(parameters_file)

    else:
        raise FileNotFoundError(f'ERROR: The parameter json file specified "{args.parameters}" is not valid.')


    if args.options: 
        if os.path.isfile(args.options) and '.json' in args.options:

            options_file_name = args.options
            print(f'The option file "{args.options}" was specified.') 

            with open(options_file_name) as options_file:
                options = json.load(options_file)

        else:
            raise FileNotFoundError(f'ERROR: The option json file specified "{args.options}" is not valid.')

    else:
        options = None
        print('No options file was specified.\nNOTE: The user will have to select the options via a cmdline interface.')


    print_postprocessing_setup(folder_to_read, options, parameters)

    post_processor = ControlPostProcessor(
        folder     = folder_to_read, 
        options    = options, 
        parameters = parameters
    )

    post_processor.post_process()


if __name__ == '__main__':
    main()
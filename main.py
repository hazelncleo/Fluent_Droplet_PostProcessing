from post_processor import PostProcessor
import sys
import os
import json

# TODO, formalise this part of the script
parameters = {
    'frequency'           : 1.63e6,
    'amplitude'           : 1e-6,
    'n_cycles'            : 60,
    'n_levels_refinement' : 4,
    'channel_width'       : 50,
    'grid_size'           : 500
}


if __name__ == '__main__':
    
    # Check for folder specification TODO: Clean this up to be bug free
    if (len(sys.argv) - 1) and os.path.isdir(sys.argv[1]):
        
        cmd_line_specified_folder = os.path.abspath(sys.argv[1])

        # Check if options .json provided
        if (len(sys.argv) - 2) and os.path.isfile(sys.argv[2]) and sys.argv[2].endswith('.json'):

            with open(sys.argv[2]) as options_file:
                options = json.load(options_file)

            print('Executing with: \nParameters = {}\nFolder = {}\nOptions = {}'.format(parameters, cmd_line_specified_folder, options))
            post_processor = PostProcessor(parameters = parameters, folder = cmd_line_specified_folder, options = options)

        else:
            print('Executing with: \nParameters = {}\nFolder = {}'.format(parameters, cmd_line_specified_folder))
            post_processor = PostProcessor(parameters = parameters, folder = cmd_line_specified_folder)

    
        
    else:
        print('Executing with: \nParameters = {}'.format(parameters))
        post_processor = PostProcessor(parameters = parameters)
    
    post_processor.post_process()
    
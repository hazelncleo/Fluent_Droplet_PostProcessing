from post_processor import PostProcessor
import sys
import os

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
    
    if (len(sys.argv) - 1) and os.path.isdir(sys.argv[1]):
        
        cmd_line_specified_folder = os.path.abspath(sys.argv[1])
        post_processor = PostProcessor(parameters = parameters, folder = cmd_line_specified_folder)
        
    else:
        post_processor = PostProcessor(parameters = parameters)
    
    post_processor.post_process()
    
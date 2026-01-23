from post_processor import PostProcessor


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
    post_processor = PostProcessor(parameters)

    post_processor.post_process()
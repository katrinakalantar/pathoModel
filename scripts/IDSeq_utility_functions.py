def parse_params(config_file):
    configs = open(config_file, 'r').readlines()
    parameters = dict([(i.strip().split(':')[0],i.strip().split(':')[1].strip()) for i in configs])
    mp = parameters['model_parameters'].split(',')
    gp = parameters['plot_x']
    meta_class = parameters['label']
    positive_class = parameters['positive_class']
    negative_class = parameters['negative_class']
    unknown_class = parameters['unknown_class'].split(',')
    training_sample_names = parameters['training_set'].split(',')
    
    output_parameters = {
        'model_parameters':mp,
        'graph_axis':gp,
        'label':meta_class,
        'positive_class':positive_class,
        'negative_class':negative_class,
        'unknown_class':unknown_class,
        'training_set':training_sample_names
    }
    
    return(output_parameters)
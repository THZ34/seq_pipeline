def shell(samples, template, **kwargs):
    """To paraphrase the shell template

    :param samples: samples name
    :param template: shell template
    :param kwargs: other parameters, base on the template
    :return: list of command
    """

    if template.__name__ in ['mageck_count', 'bwa']:
        command_list = template(samples, kwargs)
    else:
        command_list = []
        for sample in samples:
            if type(sample) == str:
                command = template(sample, kwargs)
                command_list.append(command)
            elif type(sample) == tuple:
                for rep in sample:
                    command = template(rep, kwargs)
                    command_list.append(command)
    if 'conda' in kwargs:
        command_list.insert(0, 'conda activate ' + kwargs['conda'])
        command_list.append('conda deactivate')
    return command_list

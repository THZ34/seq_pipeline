def shell(samples, template, **kwargs):
    """"""

    if template.__name__ == 'mageck_count':
        command_list = template(samples,kwargs)
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
    return command_list



# check config files only
rule check:
    input:
        config["sample_info"]                  # do nothing - this file should exist

# print out the configuration
rule showconf:
    input:
        config["sample_info"]
    run:
        import yaml
        print('# full aggregated configuration:')
        print(yaml.dump(config).strip())
        print('# END')

###

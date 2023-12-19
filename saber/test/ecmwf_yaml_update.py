#!/usr/bin/env python3
# This script creates equivalent yaml files for a use with OOPS-ECMWF

import argparse
import datetime
import os
import yaml
import copy
import sys
import collections.abc

# Correct date formatting
def correct_date(d):
    if isinstance(d, dict):
        # Loop over dictionary items
        for k,v in d.items():
            d[k] = correct_date(v)
    elif isinstance(d, list):
        # Loop over list items
        i = 0
        for v in d:
            d[i] = correct_date(v)
            i += 1

    if isinstance(d, datetime.datetime):
        # Replace with string
        d = d.strftime("%Y-%m-%dT%H:%M:%SZ")
    return d

# Replace key name recursively
def rename_key(d, key, repl):
    if isinstance(d, dict):
        # Loop over dictionary items
        newd = {}
        for k,v in d.items():
            if k == key:
                newd[repl] = v
            else:
                newd[k] = rename_key(v, key, repl)
    elif isinstance(d, list):
        # Loop over list items
        newd = []
        for v in d:
            newd.append(rename_key(v, key, repl))
    else:
        newd = d
    return newd

# Find pattern in value and replace it
def find_and_replace(d, value, repl):
    if isinstance(d, dict):
        # Loop over dictionary items
        newd = {}
        for k,v in d.items():
             newd[k] = find_and_replace(v, value, repl)
    elif isinstance(d, list):
        # Loop over list items
        newd = []
        for v in d:
            newd.append(find_and_replace(v, value, repl))
    elif isinstance(d, str):
        newd = d.replace(value, repl)
    else:
        newd = d
    return newd

# Expand ensemble template
def expand_ensemble_template(d):
    if isinstance(d, dict):
        # Loop over dictionary items
        newd = {}
        for k,v in d.items():
            if k == "members from template":
                template = v["template"]
                pattern = v["pattern"]
                members = v["nmembers"]
                if "zero padding" in v:
                    zpad = v["zero padding"]
                else:
                    zpad = 0
                if "start" in v:
                    start = v["start"]
                else:
                    start = 1
                if "except" in v:
                    excpt = v["except"]
                else:
                    excpt = []
                is4D = ("states" in template)
                if is4D:
                    n4D = len(template["states"])
                    newd["members"] = []
                else:
                    n4D = 1
                i4D = 0
                while i4D < n4D:
                    ensemble = []
                    index = start
                    member = 0
                    while member < members:
                        if not index in excpt:
                            if is4D:
                                state = copy.deepcopy(template["states"][i4D])
                            else:
                                state = copy.deepcopy(template)
                            state = find_and_replace(state, pattern, str(index).zfill(zpad))
                            ensemble.append(state)
                            member += 1
                        index += 1
                    if is4D:
                        newd["members"].append(ensemble)
                    else:
                        newd["members"] = ensemble
                    i4D += 1
            else:
                newd[k] = expand_ensemble_template(v)
    elif isinstance(d, list):
        # Loop over list items
        newd = []
        for v in d:
            newd.append(expand_ensemble_template(v))
    else:
        newd = d
    return newd

# Add ensemble variables
def add_ensemble_variables(config):
    if "ensemble" in config:
        covariance = config["ensemble"]
        ensemble = covariance["members"]
        is4D = isinstance(ensemble[0], collections.abc.Sequence)
        covariance.pop("members")
        if is4D:
            n4D = len(ensemble)
            config["ensemble"] = []
            for i4D in range(0, n4D):
                conf3D = {}
                conf3D["state"] = ensemble[i4D]
                conf3D["variables"] = ensemble[0][0]["variables"]
                conf3D["members"] = len(ensemble[0])
                config["ensemble"].append(conf3D)
        else:
            covariance["state"] = ensemble
            covariance["date"] = date
            covariance["variables"] = ensemble[0]["variables"]
            covariance["members"] = len(ensemble)

    if "ensemble pert" in config:
        covariance = config["ensemble pert"]
        ensemble = covariance["members"]
        covariance["state"] = ensemble
        covariance["date"] = date
        covariance["variables"] = variables
        covariance.pop("members")
        covariance["members"] = len(ensemble)

    if "dual resolution calibration" in config:
        dualConfig = config["dual resolution calibration"]

        if "ensemble" in dualConfig:
            covariance = dualConfig["ensemble"]
            ensemble = covariance["members"]
            covariance["state"] = ensemble
            covariance["date"] = date
            covariance["variables"] = variables
            covariance.pop("members")
            covariance["members"] = len(ensemble)

        if "ensemble pert" in dualConfig:
            covariance = dualConfig["ensemble pert"]
            ensemble = covariance["members"]
            covariance["state"] = ensemble
            covariance["date"] = date
            covariance["variables"] = variables
            covariance.pop("members")
            covariance["members"] = len(ensemble)

    if "members" in config:
        ensemble = config["members"]
        config["ensemble"] = []
        is4D = isinstance(ensemble[0], collections.abc.Sequence)
        config.pop("members")
        if is4D:
            n4D = len(ensemble)
            for i4D in range(0, n4D):
                conf3D = {}
                conf3D["state"] = ensemble[i4D]
                conf3D["variables"] = variables
                conf3D["members"] = len(ensemble[0])
                config["ensemble"].append(conf3D)
        else:
            conf3D = {}
            conf3D["state"] = ensemble
            conf3D["variables"] = variables
            conf3D["members"] = len(ensemble)
            config["ensemble"].append(conf3D)

    return config

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("inputYaml", help="Yaml file name")
parser.add_argument("outputYaml", help="Yaml file name")
args = parser.parse_args()
#print("Processing " + args.inputYaml + " into " + args.outputYaml)

# Exit for 4D files (not ready yet)
if "_4d" in args.inputYaml:
    exit()

# Read yaml file
print("--  - Update yaml: " + args.inputYaml)
with open(args.inputYaml, "r") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Correct date formatting recursively
correct_date(config)

# Rename "geometry" into "resolution"
if "geometry" in config:
    geometry = config["geometry"]
    config["resolution"] = geometry
    config.pop("geometry")

# Rename "state variables" into "variables"
config = rename_key(config, "state variables", "variables")

# Expand ensemble template
config = expand_ensemble_template(config)

# Background
if "background" in config:
    background = {}
    if "states" in config["background"]:
      background["state"] = config["background"]["states"]
    else:
      background["state"] = [config["background"]]
    config["Background"] = background
    config.pop("background")
    date = config["Background"]["state"][0]["date"]
    variables = config["Background"]["state"][0]["variables"]

# Background error
if "background error" in config:
    covariance = config["background error"]
    covariance["covariance"] = covariance["covariance model"]
    covariance.pop("covariance model")
    config["Covariance"] = covariance
    config.pop("background error")

    # Ensemble at covariance root
    config["Covariance"] = add_ensemble_variables(config["Covariance"])

    # Ensemble at SABER hybrid components root
    if config["Covariance"]["covariance"] == "SABER":
        if config["Covariance"]["saber central block"]["saber block name"] == "Hybrid":
            for component in config["Covariance"]["saber central block"]["components"]:
                component["covariance"] = add_ensemble_variables(component["covariance"])

    # OOPS hybrid case
    if config["Covariance"]["covariance"] == "hybrid":
        # Check whether all components are hybrid
        components = config["Covariance"]["components"]
        allSaber = True
        if len(components) == 2:
            if components[1]["covariance"]["covariance model"] == "ensemble":
                allSaber = False

        if allSaber:
            # Switch to SABER hybrid block (OOPS-ECMWF hybrid not flexible enough yet)
            config["Covariance"]["covariance"] = "SABER"
            centralBlock = {}
            centralBlock["saber block name"] = "Hybrid"
            centralBlock["components"] = components
            config["Covariance"]["saber central block"] = centralBlock
            for component in components:
                component["covariance"] = add_ensemble_variables(component["covariance"])
            adjTest = False
            sqrtTest = False
            for component in components:
                if "adjoint test" in component["covariance"]:
                    adjTest = adjTest or component["covariance"]["adjoint test"]
                    component["covariance"].pop("adjoint test")
                if "square-root test" in component["covariance"]:
                    sqrtTest = sqrtTest or component["covariance"]["square-root test"]
                    component["covariance"].pop("square-root test")
                if "square-root tolerance" in component["covariance"]:
                    component["covariance"]["saber central block"]["square-root tolerance"] = component["covariance"]["square-root tolerance"]
                    component["covariance"].pop("square-root tolerance")
            config["Covariance"]["adjoint test"] = adjTest
            config["Covariance"]["square-root test"] = sqrtTest

        else:
            # Update static_covariance
            covariance = components[0]["covariance"]
            static_weight = components[0]["weight"]["value"]
            if covariance["covariance model"] == "SABER":
                config["Covariance"]["static_covariance"] = covariance
            covariance = add_ensemble_variables(covariance)
            config["Covariance"]["static_weight"] = {}
            config["Covariance"]["static_weight"]["matrix"] = "hybrid_weight"
            config["Covariance"]["static_weight"]["weight"] = static_weight            

            # Update ensemble
            covariance = components[1]["covariance"]
            ensemble_weight = components[1]["weight"]["value"]
            if covariance["covariance model"] == "ensemble":
                if "localization" in covariance:
                    covariance["localization"]["localization"] = covariance["localization"]["localization method"]
                    covariance["localization"].pop("localization method")
                    covariance["localization"]["variables"] = variables
                covariance = add_ensemble_variables(covariance)
                config["Covariance"]["ensemble_covariance"] = covariance
                config["Covariance"]["ensemble_weight"] = {}
                config["Covariance"]["ensemble_weight"]["matrix"] = "hybrid_weight"
                config["Covariance"]["ensemble_weight"]["weight"] = ensemble_weight

            # Update covariance model
            config["Covariance"]["static_covariance"]["covariance"] = config["Covariance"]["static_covariance"]["covariance model"]
            config["Covariance"]["static_covariance"].pop("covariance model")
            config["Covariance"]["ensemble_covariance"]["covariance"] = config["Covariance"]["ensemble_covariance"]["covariance model"]
            config["Covariance"]["ensemble_covariance"].pop("covariance model")

        # Remove components
        config["Covariance"].pop("components")

    # Remove linear variable change
    if "linear variable change" in config["Covariance"]:
        config["Covariance"].pop("linear variable change")

# Ensemble
config = add_ensemble_variables(config)

# Write yaml file
with open(args.outputYaml, "w") as file:
    output = yaml.dump(config, file, sort_keys=False)

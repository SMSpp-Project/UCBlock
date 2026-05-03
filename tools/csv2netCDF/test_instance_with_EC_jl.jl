## Test the instances with EnergyCommunity.jl
#
# This file aims to obtain the results of the instances with EnergyCommunity.jl.
#
# To run this file you can run in the terminal:
# julia test_instance_with_EC_jl.jl
# the results will be printed in the terminal.
#
# To customize the instance, you can change the configuration file (variable fconfig).
#
# Each instance is defined by:
# - configuration file: configuration file of the instance, which contains the parameters of the model. 
#   The configuration file is a YAML file, e.g. `energy_community_model.yaml`
# - data file(s): data file of the instance, which contains the data of the model. 
#   The data file(s) are CSV files, e.g. `input_resources.csv` as defined in the 
#   configuration file under general->optional_datasets
#
# For more details on the configuration files, please see the documentation of EnergyCommunity.jl:
# https://spsunipi.github.io/EnergyCommunity.jl/dev/configuration/configuration/
#
# This file is structured as follows:
# 1. Imports the necessary packages
# 2. Defines the path to the configuration file
# 3. Create the model and solve it
# 4. Print the results

## 1. Imports the necessary packages

# Package manager to setup the environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()  # optional; comment after first execution

using EnergyCommunity
using Gurobi  # Commercial solver; if you don't have a license you can use HiGHS
using HiGHS
using JuMP

# default solver
optimizer = Gurobi.Optimizer  # if you have a Gurobi license, otherwise use HiGHS
# optimizer = HiGHS.Optimizer  # Uncomment this line to use HiGHS instead of Gurobi

## 2. Defines the path to the configuration file

# Load the data
fconfig = "energy_community_model_new.yml"

## 3. Create the model and solve it

# Create the model
model = ModelEC(fconfig, EnergyCommunity.GroupCO(), optimizer)

# build the model
build_model!(model)

# Solve the model
optimize!(model)

## 4. Print the results

println("Optimal value: ", objective_value(model))
println("Optimal installed capacity by user: ", value.(model.results[:x_us]))
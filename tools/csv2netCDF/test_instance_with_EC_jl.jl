## Test the instances with EnergyCommunity.jl
#
# This file aims to obtain the results of the instances with EnergyCommunity.jl.
# If a stochasic instance is requested, it will be solved with the stochastic version of EnergyCommunity.jl, otherwise it will be solved with the deterministic version of EnergyCommunity.jl.
# A stochastic instance is identified by the presence of "_sto.yml" in the name of the configuration file, otherwise it is considered deterministic.
#
# To run this file you can run in the terminal:
# julia test_instance_with_EC_jl.jl {file_name [optional]}
# the results will be printed in the terminal.
#
# {file_name} is the optional name of the configuration file, which is a YAML file, e.g. `energy_community_model.yaml`.
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
fconfig = "energy_community_model_new.yml"  # default value
if length(ARGS) > 0
    fconfig = string(ARGS[1])
end
println("Using configuration file: ", fconfig)

# define if network is stochastic: if the configuration file contains "_sto.yml" it is considered stochastic, otherwise it is deterministic

is_stochastic = "_sto.yml" in fconfig

# Ensure environment is set up with the correct version of EnergyCommunity.jl
if is_stochastic in fconfig
    Pkg.add(url="https://github.com/SPSUnipi/EnergyCommunity.jl", rev="stochastic")
else
    Pkg.add(url="https://github.com/SPSUnipi/EnergyCommunity.jl", rev="main")
end

## 3. Create the model and solve it

obj_value = nothing
optimal_design = nothing

if is_stochastic
    println("The model is stochastic.")

    # TBD
else
    println("The model is deterministic.")

    # Create the model
    model = ModelEC(fconfig, EnergyCommunity.GroupCO(), optimizer)

    # build the model
    build_model!(model)

    # Solve the model
    optimize!(model)

    obj_value = objective_value(model)
    optimal_design = value.(model.results[:x_us])
end


## 4. Print the results

println("Optimal value: ", obj_value)
println("Optimal installed capacity by user: ", optimal_design)
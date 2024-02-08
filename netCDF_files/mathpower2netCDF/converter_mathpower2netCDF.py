# author : Quentin Jacquet
import sys, os, subprocess, argparse
import netCDF4
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir,'..'))
from format_matpower2netCDF import *

class ConverterMathpower2netCDF:
    def __init__(self, filepath, mode = 'DC'):
        """
        Parse the .m file and fill self.attrs such that:
        self.attrs = {"mpc.key": [line for line in range(nb_of_key)]}
        """
        self.str_file = None
        self.mode = mode
        self.attrs = {}

        with open(filepath,'r') as f:
            self.str_file = f.read().replace(" ", "")

        # name of the instance on first line
        self.attrs["mpc.name"] = self.str_file[self.str_file.find("=")+1:self.str_file.find("\n")]
        
        # split with end of line
        self.str_file = self.str_file.split(";") # matlab file: lines end with ';'
        
        # find start of line
        for i,l in enumerate(self.str_file):
            start_l = l.find("mpc.")
            if start_l >= 0:
                self.str_file[i] = l[start_l:].split("\t")
            else:
                self.str_file[i] = l.split("\t")
        
        # fill data (either single line, or tabular of data)
        i = 0
        while i < len(self.str_file):
            l  = self.str_file[i]
            key = l[0][:l[0].find("=")]

            # if data given by a tabular
            if "[" in self.str_file[i][0]:
                tab = []
                j = i
                while "]" not in self.str_file[j][-1]: # assumption: the tab ends with a line "];"
                    if "\n" in self.str_file[j][0]:
                        self.str_file[j] = [float(k) for k in self.str_file[j][1:]]
                    else:
                        self.str_file[j] = [float(k) for k in self.str_file[j]]
                    tab.append(self.str_file[j])
                    j += 1
                i = j + 1
                self.attrs[key] = tab

            # if data is a just a value
            else:
                val = l[0][l[0].find("=")+1:]
                try:
                    self.attrs[key] = float(val)
                except:
                    self.attrs[key] = val
                i += 1
        print("INFO: Input file '{0}' has been correcty read".format(filepath))


    def create_nc_file(self, filepath, txt_output = False):
        """
        Write the .nc file corresponding the instance previously loaded.
        If 'txt_output' is True, the exe 'ncdump.exe' has to be defined in the env variables.
        """
        rootgrp = netCDF4.Dataset("{0}.nc".format(filepath), "w", format="NETCDF4")
        rootgrp.setncatts({'SMS++_file_type':1}) # global attribute
        
        maingrp = rootgrp.createGroup("Block_0") # main block
        dimensions = {  "mpc.gen":      len(self.attrs["mpc.gen"]),
                        "mpc.bus":      len(self.attrs["mpc.bus"]),
                        "mpc.branch":   len(self.attrs["mpc.branch"])
        }
        for mpc_lab, d in dimensions.items():
            maingrp.createDimension(dim_labels[mpc_lab], d)

        # subgroups by generator
        subgroups = []
        for i in range(dimensions["mpc.gen"]):
            subgroups.append(maingrp.createGroup("UnitBlock_{0}".format(i)))
            subgroups[i].setncatts({'type':"ThermalUnitBlock"})

        # create variables and fill data
        for mpc_lab in ["mpc.gen","mpc.bus","mpc.branch"]:
            lab_tab = var_labels[mpc_lab]
            for idx,t in enumerate(lab_tab):
                k,v,*g = t # k = label of the netCDF list, v is the type (or None), len(g) > 0 means that it is in generator subgroup
                
                # for main group
                if v is not None and len(g) == 0:
                    var = maingrp.createVariable(k,v,dim_labels[mpc_lab])
                    if "int" in v:      var[:] = [int(l[idx]) for l in self.attrs[mpc_lab]]
                    if "double" in v:   var[:] = [float(l[idx]) for l in self.attrs[mpc_lab]]
                # for subgroups
                if v is not None and len(g) > 0:
                    for i,l in enumerate(self.attrs[mpc_lab]):
                        var = subgroups[i].createVariable(k,v,chunksizes=(1,))
                        if "int" in v:      var[:] = int(l[idx])
                        if "double" in v:   var[:] = float(l[idx])

        # add costs in each subgroup ThermalUnitBlock
        mpc_lab = "mpc.gencost"
        if mpc_lab in self.attrs:
            lab_tab = var_labels[mpc_lab]
            idx = 0
            gen_dim = dimensions["mpc.gen"]
            while idx < len(lab_tab):
                k,v,*g = lab_tab[idx] # k = label of the netCDF list, v is the type (or None), len(g) > 0 means that it is in generator subgroup
                if v is not None and len(g) > 0:
                    for i,l in enumerate(self.attrs[mpc_lab]):
                        if i < gen_dim:
                            var = subgroups[i].createVariable(k,v,chunksizes=(1,))
                        elif i < 2*gen_dim: # for reactive coeffs (optional)
                            var = subgroups[i%gen_dim].createVariable("Reactive{0}".format(k),v,chunksizes=(1,))
                        if "int" in v:      var[:] = int(l[idx])
                        if "double" in v:   var[:] = float(l[idx])
                idx += 1

            # add costs coefficients
            for i,l in enumerate(self.attrs[mpc_lab]):
                if i < gen_dim:
                    nb_coeff_cost = int(l[idx])
                    subgroups[i].createDimension("NumberCostCoeffs", nb_coeff_cost)
                    var = subgroups[i].createVariable("PowerCostCoeffs", "double", "NumberCostCoeffs")
                    var[:] = [float(l[idx + k + 1]) for k in range(nb_coeff_cost)]
                elif i < 2*gen_dim: # for reactive coeffs (optional)
                    nb_coeff_cost = int(l[idx])
                    subgroups[i%gen_dim].createDimension("NumberReactiveCostCoeffs", nb_coeff_cost)
                    var = subgroups[i%gen_dim].createVariable(  "ReactivePowerCostCoeffs", "double", "NumberReactiveCostCoeffs")
                    var[:] = [float(l[idx + k + 1]) for k in range(nb_coeff_cost)]

        rootgrp.close()
        print("INFO: File '{0}.nc' written".format(filepath))

        if txt_output:
            with open("{0}.txt".format(filepath), "w") as txt_file:
                p = subprocess.Popen(["ncdump.exe", "{0}.nc".format(filepath)], stdout=txt_file, stderr=subprocess.PIPE)
            

if __name__.endswith("__main__"):
    parser = argparse.ArgumentParser(   prog='ConverterMathpower2netCDF',
                                        description='Converter Mathpower -> netCDF')
    parser.add_argument('folder', metavar = "<input_folder>", type=str,
                        help = 'input folder path')
    parser.add_argument('-t', '--type', choices = ['AC', 'DC'], default = 'DC',
                        help = 'type of instance')
    parser.add_argument('-f', '--format', choices = ['txt', 'nconly'], default = 'nconly',
                        help = "format of the output")

    args = parser.parse_args()

    for root, dirs, files in os.walk(".", topdown=False):
        for i, name in enumerate(files):
            filepath = os.path.join(root, name)
            if filepath.endswith(".m"):
                try:
                    output_filename = "{0}".format(filepath[:-2], args.type)
                    converter = ConverterMathpower2netCDF(filepath, mode = args.type)
                    converter.create_nc_file(output_filename, txt_output = (args.format == "txt"))
                except:
                    print("Somethings goes wrong with '{0}'".format(filepath))
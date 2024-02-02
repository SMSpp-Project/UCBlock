# author : Quentin Jacquet
import sys,os
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir,'..'))
from format_matpower2netCDF import *

class ConverterMathpower2netCDF:
    def __init__(self, file_name):
        self.str_file = None
        self.attrs = {}

        with open(file_name,'r') as f:
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
            else:
                val = l[0][l[0].find("=")+1:]
                try:
                    self.attrs[key] = float(val)
                except:
                    self.attrs[key] = val
                i += 1


    def create_nc_file(self, file_name):
        block = "group: Block_0 {"
        
        # ===== DIMENSIONS
        dimensions = {
            "NumberElectricalGenerators": len(self.attrs["mpc.gen"]),
            "NumberNodes": len(self.attrs["mpc.bus"]),
            "NumberLines": len(self.attrs["mpc.branch"])
        }
        block += "\ndimensions:"
        for k,v in dimensions.items():
            block += "\n\t{0} = {1} ;".format(k,v)

        # ===== VARIABLES
        block += "\n\nvariables:"

        # Nodes/Bus
        for k,v in labels_bus:
            if v is not None:
                block += "\n\t{0} {1}(NumberNodes) ;".format(v,k)
        # Generators
        for k,v in labels_gen:
            if v is not None:
                block += "\n\t{0} {1}(NumberElectricalGenerators) ;".format(v,k) 
        # Branches/Lines
        for k,v in labels_branch:
            if v is not None:
                block += "\n\t{0} {1}(NumberLines) ;".format(v,k) 

        block += '\n\n\t// group attributes:'
        block += '\n\t\t:id = "0" ;'
        block += '\n\t\t:type = "UCBlock" ;'

        # ===== DATA
        block += "\n\ndata:"
        
        # info buses
        for idx,t in enumerate(labels_bus):
            k,v = t
            if v is not None:
                block += "\n\t{0} = ".format(k)
                for i,l in enumerate(self.attrs["mpc.bus"]):
                    block += "{0},".format(int(l[idx]))
                block = block[:-1] + ";"
        # info generator
        for idx,t in enumerate(labels_gen):
            k,v = t
            if v is not None:
                block += "\n\t{0} = ".format(k)
                for i,l in enumerate(self.attrs["mpc.gen"]):
                    block += "{0},".format(int(l[idx]))
                block = block[:-1] + ";"
        # info branch
        for idx,t in enumerate(labels_branch):
            k,v = t
            if v is not None:
                block += "\n\t{0} = ".format(k)
                for i,l in enumerate(self.attrs["mpc.branch"]):
                    block += "{0},".format(int(l[idx]))
                block = block[:-1] + ";"
        # fermeture du block
        block += "\n}"


        # ===== WRITING FILE        
        with open(file_name, "w") as f:
            f.write("netcdf " + self.attrs["mpc.name"] + " {\n:SMS++_file_type = 1 ;\n" + block + "\n}")
        

if __name__.endswith("__main__"):
    if len(sys.argv) == 1:
        print("You must provide a matlab file name (.m)")
        sys.exit()
    file_name = sys.argv[1]
    if len(sys.argv) > 2:   output_file_name = sys.argv[2]
    else:                   output_file_name = file_name[:file_name.rfind(".m")]+".txt"
    converter = ConverterMathpower2netCDF(file_name)
    converter.create_nc_file(output_file_name)
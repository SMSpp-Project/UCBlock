# author : Quentin Jacquet
import sys,os,argparse
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir,'..'))
from format_matpower2netCDF import *

class ConverterMathpower2netCDF:
    def __init__(self, file_name, mode = 'DC'):
        self.str_file = None
        self.mode = mode
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
            "mpc.gen":      len(self.attrs["mpc.gen"]),
            "mpc.bus":      len(self.attrs["mpc.bus"]),
            "mpc.branch":   len(self.attrs["mpc.branch"]),
            "mpc.gencost":  int(self.attrs["mpc.gencost"][0][3]) 
            # assumption: same number of coeffs for each generator cost function
        }
        block += "\ndimensions:"
        for k,v in dimensions.items():
            block += "\n\t{0} = {1} ;".format(dim_labels[k],v)

        # ===== VARIABLES
        block += "\n\nvariables:"
        for l1,l2 in dim_labels.items():
            for k,v in labels[l1]:
                if v is not None:
                    block += "\n\t{0} {1}({2}) ;".format(v,k,l2) 
        block += "\n\tdouble PowerCostCoeffs({0},{1}) ;".format(dim_labels["mpc.gen"],dim_labels["mpc.gencost"])

        # Reference node
        block += "\n\tuint ReferenceNode ;"

        block += '\n\n\t// group attributes:'
        block += '\n\t\t:id = "0" ;'
        block += '\n\t\t:type = "{0}NetworkBlock" ;'.format(self.mode)

        # ===== DATA
        block += "\n\ndata:"
        
        # info buses, generators and branches
        for l_name, l_tab in labels.items():
            for idx,t in enumerate(l_tab):
                k,v = t # k = label of the netCDF list, v is the type (or None)
                if v is not None:
                    block += "\n\t{0} = ".format(k)
                    for i,l in enumerate(self.attrs[l_name]):
                        if "int" in v:
                            block += "{0},".format(int(l[idx]))
                        if "double" in v:
                            block += "{0},".format(float(l[idx]))
                    block = block[:-1] + ";"

        # Cost coeffs for generator (assumption: no reactive power cost, only active power cost)
        block += "\n\tPowerCostCoeffs = "
        for i,l in enumerate(self.attrs["mpc.gencost"]):
            for idx in range(4,4+dimensions["mpc.gencost"]):
                block += "{0},".format(float(l[idx]))
        block = block[:-1] + ";"

        # for reference node
        for i,l in enumerate(self.attrs["mpc.bus"]):
            if l[1] == 3:
                block += "\n\tReferenceNode = {0} ;".format(int(l[0]))
                break

        # fermeture du block
        block += "\n}"


        # ===== WRITING FILE        
        with open(file_name, "w") as f:
            f.write("netcdf " + self.attrs["mpc.name"] + " {\n:SMS++_file_type = 1 ;\n" + block + "\n}")
        

if __name__.endswith("__main__"):
    parser = argparse.ArgumentParser(
                    prog='ConverterMathpower2netCDF',
                    description='Converter Mathpower -> netCDF',
                    epilog='')
    parser.add_argument('filename', metavar = "<input>.m", type=str,
                        help = 'input file path')
    parser.add_argument('-o', '--output', metavar = "<output>.txt", type=str,
                        help = 'output file path (default: <input>.txt)')
    parser.add_argument('-t', '--type', choices = ['AC', 'DC'], default = 'DC',
                        help = 'type of instance')

    args = parser.parse_args()
    output_filename = args.output
    if output_filename is None:
        output_filename = "{0}_{1}.txt".format(args.filename[:args.filename.rfind(".m")], args.type)
    converter = ConverterMathpower2netCDF(args.filename, mode = args.type)
    converter.create_nc_file(output_filename)
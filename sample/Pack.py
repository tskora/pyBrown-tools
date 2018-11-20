from ReadInput import InputData
from Grid import pack_molecules
from WriteStructure import write_structure
from ParseInput import parse_input_filename

input_filename = parse_input_filename()
i = InputData(input_filename, [])
print(i.input_data)
coords = pack_molecules(i)
print(coords)
write_structure(i, coords)
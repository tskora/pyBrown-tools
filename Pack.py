from ReadInput import InputData
from Grid import pack_molecules

i = InputData("input.json", [])
print(i.input_data)
pop_grid = pack_molecules(i)
print(pop_grid)
def read_energies(input_data):

	input_enr_filenames = [ input_data["input_enr_template"] + str(i) + '.enr'
							for i in range( *input_data["input_enr_range"] ) ]

	number_of_enr_files = len( input_enr_filenames )

	numbers_of_timeframes = [ _file_length(input_enr_filenames[i]) - 1
							  for i in range( number_of_enr_files ) ]

	number_of_timeframes = min(numbers_of_timeframes)

	energies = []
	times = []

	for i, input_enr_filename in enumerate(input_enr_filenames):

		if input_data["verbose"]: timestamp('enr file {}/{}'.format( i + 1, len(input_enr_filenames) ))

		with open(input_enr_filename, 'r') as input_enr_file:

			_read_energies_from_enr_file(input_enr_file, energies, times)

	if input_data["verbose"]: timestamp('enrs read')

	for i in range(number_of_enr_files):

		energies[i] = energies[i][:number_of_timeframes]

	times = times[:number_of_timeframes]

	return np.array(energies), np.array(times, float)




def compute_menergies(energies):

	menergies = np.zeros( len( energies[0] ), float )

	for energies_one_file in energies:

		menergies += energies_one_file

	return menergies / len( energies )




def plot_menergies(input_data, times, menergies):

	plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'E [$\frac{kcal}{mol}$]')

	plt.plot(times / 1000000, menergies, '-')

	plt.savefig(input_data["input_enr_template"] + 'enr.jpg', dpi = 100)

	plt.close()




def _read_energies_from_enr_file(enr_file, energies, times):

	energies_from_current_file = []

	if len(times) == 0:
		add_times = True
	else:
		add_times = False
		len_i = len(times)

	enr_file.readline()

	for i, line in enumerate(enr_file):

		if add_times:
			times.append( float( line.split()[0] ) )

		else:
			if i >= len_i:
				break
			else:
				assert( float( line.split()[0] ) == times[i] )

		energies_from_current_file.append( float( line.split()[1] ) )

	energies.append( np.array( energies_from_current_file, float ) )


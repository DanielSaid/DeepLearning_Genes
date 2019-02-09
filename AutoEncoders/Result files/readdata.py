import os
import scipy.io as sio


def main():
	# path = os.path.join("Rat_vitro_human_vitro", "Random_variance")
	path = os.path.join("Rat_vitro_rat_vivo", "Random_variance")

	file = os.path.join(path, 'ae_random_2.mat')

	data = sio.loadmat(file)
	input_genes = data['genes_input']
	output_genes = data['genes']

	print("input genes:")
	for g in input_genes:
		print(g.strip())

	print("output genes:")
	for g in output_genes:
		print(g.strip())


if __name__ == '__main__':
	main()
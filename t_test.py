from scipy.stats import ttest_rel
import numpy as np
import os
import scipy.io as sio
import pandas as pd


def t_test(errors):
	p_matrix = np.zeros((errors.shape[0], errors.shape[0]))

	for i in range(errors.shape[0]):
		for j in range(errors.shape[0]):
			t_score, p_val = ttest_rel(errors[i], errors[j])
			p_matrix[i, j] = p_val

	return p_matrix


def read_uda():
	files = ['SMALL_50_SourceValidationMAEs.npy', 'CH', 'GTX_CAR_SourceValidationMAEs.npy', 'NAFLD', 'R1', 'R2']

	vivo_errs = []

	for filename in files:
		if filename[-3:] == 'npy':
			vivo_errs.append(np.load(os.path.join('UDA T-test', filename)))
		else:
			vivo_errs.append(np.zeros((720, )))

	return vivo_errs


def read_alex_ae():
	root = os.path.join('Alex T-test')
	files = ['ST', 'CH', 'GTX_CAR', 'NAFLD', 'R1', 'R2']

	vitro_errs = []
	vivo_errs = []

	for filename in files:
		vitro_data = sio.loadmat(os.path.join(root, 'vitro_alex_' + filename))
		vitro_errs.append(np.abs(vitro_data['Y_valid'] - vitro_data['recon_valid']).mean(axis=1))

		vivo_data = sio.loadmat(os.path.join(root, 'vivo_alex_' + filename))
		vivo_errs.append(np.abs(vivo_data['Y_valid'] - vivo_data['recon_valid']).mean(axis=1))

	return vitro_errs, vivo_errs


def read_ae_and_rf():
	root = os.path.join('AutoEncoders', 'Result files')
	vitro_root = os.path.join(root, 'Rat_vitro_human_vitro')
	vivo_root = os.path.join(root, 'Rat_vitro_rat_vivo')
	files = ['50', 'cholestasis', 'gtx_car', 'nafld']
	random_files = ['random_1', 'random_2']
	vitro_random = os.path.join(vitro_root, 'Random_variance')
	vivo_random = os.path.join(vivo_root, 'Random_variance')

	rf_vitro_errs = []
	ae_vitro_errs = []

	rf_vivo_errs = []
	ae_vivo_errs = []

	# main 4 sets
	for filename in files:
		# Vitro
		rf_name = os.path.join(vitro_root, 'rf_' + filename + '.mat')
		ae_name = os.path.join(vitro_root, 'ae_' + filename + '.mat')

		rf_data = sio.loadmat(rf_name)
		rf_vitro_errs.append(np.abs(rf_data['Y_valid'] - rf_data['recon_valid']).mean(axis=1))
		ae_vitro_errs.append(np.abs(rf_data['Y_valid'] - sio.loadmat(ae_name)['recon_valid']).mean(axis=1))

		# Vivo
		rf_name = os.path.join(vivo_root, 'rf_' + filename + '.mat')
		ae_name =  os.path.join(vivo_root, 'ae_' + filename + '.mat')

		rf_data = sio.loadmat(rf_name)
		rf_vivo_errs.append(np.abs(rf_data['Y_valid'] - rf_data['recon_valid']).mean(axis=1))
		ae_vivo_errs.append(np.abs(rf_data['Y_valid'] - sio.loadmat(ae_name)['recon_valid']).mean(axis=1))

	# random sets
	for filename in random_files:
		# Vitro
		rf_name = os.path.join(vitro_random, 'rf_' + filename + '.mat')
		ae_name =  os.path.join(vitro_random, 'ae_' + filename + '.mat')

		rf_data = sio.loadmat(rf_name)
		rf_vitro_errs.append(np.abs(rf_data['Y_valid'] - rf_data['recon_valid']).mean(axis=1))
		ae_vitro_errs.append(np.abs(rf_data['Y_valid'] - sio.loadmat(ae_name)['recon_valid']).mean(axis=1))

		# Vivo
		rf_name = os.path.join(vivo_random, 'rf_' + filename + '.mat')
		ae_name =  os.path.join(vivo_random, 'ae_' + filename + '.mat')

		rf_data = sio.loadmat(rf_name)
		rf_vivo_errs.append(np.abs(rf_data['Y_valid'] - rf_data['recon_valid']).mean(axis=1))
		ae_vivo_errs.append(np.abs(rf_data['Y_valid'] - sio.loadmat(ae_name)['recon_valid']).mean(axis=1))

	return rf_vitro_errs, rf_vivo_errs, ae_vitro_errs, ae_vivo_errs


def read_bae():
	root = os.path.join('BAE T-test')
	files = ['ST', 'CH', 'GTX_CAR', 'NAFLD', 'R1', 'R2']

	vitro_errs = []
	vivo_errs = []

	# main 4 sets
	for filename in files:
		vitro_data = np.load(os.path.join(root, filename + '_vitro.npy')).item()
		vitro_errs.append(np.abs(vitro_data['Y_valid'] - vitro_data['recon_valid']).mean(axis=1))

		vivo_data = np.load(os.path.join(root, filename + '_vivo.npy')).item()
		vivo_errs.append(np.abs(vivo_data['Y_valid'] - vivo_data['recon_valid']).mean(axis=1))

	return vitro_errs, vivo_errs


def read_cnn():
	root = os.path.join('CNN T-test')
	files = ['st.mat', 'ch.mat', 'gtx_car.mat', 'nafld.mat', 'r1.mat', 'r2.mat']

	vitro_errs = []
	vivo_errs = []

	for filename in files:
		vitro_data = sio.loadmat(os.path.join(root, 'cnn_vitro_' + filename))
		vitro_errs.append(np.abs(vitro_data['Y_valid'] - vitro_data['recon_valid']).mean(axis=1))

		vivo_data = sio.loadmat(os.path.join(root, 'cnn_vivo_' + filename))
		vivo_errs.append(np.abs(vivo_data['Y_valid'] - vivo_data['recon_valid']).mean(axis=1))

	return vitro_errs, vivo_errs


def main():
	rf_vitro_errs, rf_vivo_errs, ae_vitro_errs, ae_vivo_errs = read_ae_and_rf()
	bae_vitro_errs, bae_vivo_errs = read_bae()
	alex_vitro_errs, alex_vivo_errs = read_alex_ae()
	cnn_vitro_errs, cnn_vivo_errs = read_cnn()
	uda_vivo_errs = read_uda()

	vitro_errors = np.array([ae_vitro_errs, rf_vitro_errs, bae_vitro_errs, alex_vitro_errs, cnn_vitro_errs])
	vivo_errors = np.array([ae_vivo_errs, rf_vivo_errs, bae_vivo_errs, alex_vivo_errs, cnn_vivo_errs, uda_vivo_errs])

	datasets = ['ST', 'CH', 'GTX_CAR', 'NAFLD', 'R1', 'R2']
	methods = ['AE', 'RF', 'BAE', 'ALEX', 'CNN', 'UDA']

	# Vitro
	print("VITRO P-VALUES:")
	for i in range(len(datasets)):
		print(datasets[i])
		p_matrix = t_test(vitro_errors[:, i, :])
		print(pd.DataFrame(data=p_matrix, columns=methods[:-1], index=methods[:-1]).to_string())
		print()

	# Vivo
	print("VIVO P-VALUES:")
	for i in range(len(datasets)):
		print(datasets[i])
		p_matrix = t_test(vivo_errors[:, i, :])
		print(pd.DataFrame(data=p_matrix, columns=methods, index=methods).to_string())
		print()



if __name__ == '__main__':
	main()
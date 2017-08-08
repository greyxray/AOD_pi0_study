import os
import subprocess
import json
from random import choice
from string import ascii_uppercase, digits, ascii_lowercase
import ahadd

# def __getRandomRootName(dir):
#     """ Return a random file name """
#     chars = ascii_uppercase + digits + ascii_lowercase
#     ran = []
#     for x in xrange(6):
#         ran.append(choice(chars))
#     num = self.counter
#     self.counter += 1
#     return dir + '/input_%i_' % (num) + ''.join(ran) + ".root"

# def lounchHadd(outFile, startNum, endNum, flist, n_join=20, depth = 0):
# 	if endNum - startNum + 1 <= n_join:
# 		os.system('hadd -f ' + outFile + " "  + ' '.join(flist[startNum:endNum]))
# 	else:
# 		os.system('mkdir input_depth_'+str(depth))
# 		for i in xrange(0, totalNum, n_join):
# 			outFile_new = outFile.split('.')[0] + depth + '.root'
# 		os.system('cd ..')

# 	# Loop to fill tuples
# 		nAtOnce = 20
# 		writeDir = "."
#         totalNum = len(crab_ds[key]['in_files'])
#         fileTuple = []
#         for i in xrange(0, totalNum, nAtOnce):
#             outFile = __getRandomRootName(writeDir)
#             tmpList = crab_ds[key]['in_files'][i:i + nAtOnce]
#             startNum = i + 1
#             endNum = i + len(tmpList)
#             newTuple = (outFile, startNum, endNum, totalNum, tmpList)
#             fileTuple.append(newTuple)

#         return tuple(fileTuple)


##### START OF CODE
if __name__ == '__main__':

	#full_data = 'dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/ohlushch/JetHT/CRAB3_tutorial_May2015_Data_analysis_JetHTdata_full/170516_113921/'
	#full_mc = 'dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/ohlushch/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/CRAB3_tutorial_May2015_Data_analysis_QCD_Pt_20toInf_MuEnrichedPt15/170516_113031/'
	# full_data = '/net/scratch_cms3b/hlushchenko/crab_outputs_K/JetHT/CRAB3_tutorial_May2015_Data_analysis_JetHTdata_full/170516_113921/'
	# full_mc = '/net/scratch_cms3b/hlushchenko/crab_outputs_K/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/CRAB3_tutorial_May2015_Data_analysis_QCD_Pt_20toInf_MuEnrichedPt15/170516_113031/'

	# extended binning
	#full_data = 'dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/ohlushch/JetHT/CRAB3_tutorial_May2015_Data_analysis_JetHTdata_finer_binning_2l/170518_090102/'
	#full_mc = 'dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/ohlushch/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/CRAB3_tutorial_May2015_Data_analysis_QCD_Pt_20toInf_MuEnrichedPt15_finer_binning_2/170518_090037/'

	# with bins of eta pt
	full_data = 'dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/ohlushch/JetHT/CRAB3_tutorial_May2015_Data_analysis_JetHTdata_HPS_matched/170807_223108/'
	full_mc = 'dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/ohlushch/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/CRAB3_tutorial_May2015_Data_analysis_QCD_Pt_20toInf_MuEnrichedPt15_HPS_matched/170807_222620/'
	# TODO: recursive under the dir
	crab_ds = {
		'data':
			{
			'in_path': [full_data + '0000/', full_data + '0001/', full_data + '0002/'],
			'out_name': 'out_data_full_finer_binning_new.root',
			'in_files': []
			},
		'mc':
			{
			'in_path': [full_mc + '0000/', full_mc + '0001/'],
			'out_name': 'out_mc_full_finer_binning_new.root',
			'in_files': []
			}
		}
	for key in ['mc']:#['data']:#crab_ds.keys():
		print "Key:", key
		for path in crab_ds[key]['in_path']:
			print "\tpath:", path, "..."
			p = os.popen('gfal-ls '+ path, "r")
			while 1:
				line = p.readline()
				if not line: break

				if line[-6:-1] != ".root": continue
				crab_ds[key]['in_files'].append(path+line[:-1])

				# # To remove unwanted hist from the rootfile
				# print 'rootrm '+ path+line[:-1] + ":h_Ks_v0_vz"
				# os.system('rootrm '+ path+line[:-1] + ":h_Ks_v0_vz")
			print "done:"

		if len(crab_ds[key]['in_files']) < 20:
			print "less then I thought"
			os.system('hadd -f '+ crab_ds[key]['out_name'] + " "  + ' '.join(crab_ds[key]['in_files']))
		else:
			ahadd.main(['-j', '1',
				'-t', '/net/scratch_cms3b/hlushchenko/crab_outputs_K/temp_merging_dir', #'/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_7_4_7/src/tmp_dir',
				'-V', '-f',# '-s',
				crab_ds[key]['out_name']] + crab_ds[key]['in_files']) #? didn't work with  *crab_ds[key]['in_files']

		print "\toutputfile:", crab_ds[key]['out_name']

	'''
	path = 'dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/ohlushch/JetHT/CRAB3_tutorial_May2015_Data_analysis/170515_213813/0000/'
	p = os.popen('gfal-ls '+ path,"r")

	files = []
	while 1:
		line = p.readline()
		if not line:
			break
		if len(line)>0:
			#print len(line), line[:-1]
			files.append(path+line[:-1])

	with open('onlyfiles.json', 'w') as outfile:
	    json.dump({'files': ' '.join(files)}, outfile)
	'''

	#print "-i "+ ' '.join(files)
	#higgsplot.HiggsPlotter(list_of_config_dicts=["-i "+' '.join(files), "-x h_Ks_v0_inv_m_pi", "--live", "--x-bins 100,0,1"])
	#list_of_config_dicts=plot_configs, list_of_args_strings=[args.args], n_processes=args.n_processes, n_plots=args.n_plots)

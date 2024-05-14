# -*- coding: utf-8 -*-
import os
import pickle
import numpy as np
from Tool import *
from time import gmtime, strftime
from os import path



class FunctionConfig:

    def file2config(self, path, config):

        with open(path, 'r') as f:

            for line in f.readlines():
                if line.startswith("#") or line.startswith("["):
                    continue

                p_EqualSign = line.find('=')

                if -1 == p_EqualSign:
                    pass
                else:
                    subLine = line.split(';')[0]
                    self.__soldierParseLine(line, config)

        with open("ini/ThermoScore.ini", 'r') as f1:

            for line in f1.readlines():

                p_EqualSign = line.find('=')

                if -1 == p_EqualSign:
                    pass
                else:
                    self.__soldierThermoRadarInfo(line, config)

        with open("ini/BrukerScore.ini", 'r') as f2:

            for line in f2.readlines():

                p_EqualSign = line.find('=')

                if -1 == p_EqualSign:
                    pass
                else:
                    self.__soldierBrukerRadarInfo(line, config)

        with open("ini/SCIEXScore.ini", 'r') as f2:

            for line in f2.readlines():

                p_EqualSign = line.find('=')

                if -1 == p_EqualSign:
                    pass
                else:
                    self.__soldierSCIEXRadarInfo(line, config)

    def __soldierParseLine(self, line, cfg):

        if line.startswith("#") or line.startswith("["):
            pass
        else:

            str_name = toolGetWord(line, 0, '=')
            str_value = toolGetWord(line, 1, '=').replace("\n", "")

            if "WORK_FLOW" == str_name:
                cfg.A0_TYPE_FLOW = int(str_value)

            elif "TYPE_DATA" == str_name:
                cfg.A1_TYPE_DATA = int(str_value)

            elif "PATH_RAW" == str_name:
                cfg.A0_PATH_RAW = str_value.split('|')[0]

            elif "TYPE_ANALYSIS_RESULT" == str_name:
                cfg.B3_TYPE_IDENTIFICATION_RESULT = int(str_value)

            elif "PATH_ANALYSIS_RESULT" == str_name:
                cfg.B2_PATH_IDENTIFICATION_RESULT = str_value

            elif "THRESHOLD_FDR" == str_name:
                cfg.B4_THRESHOLD_FDR = float(str_value)

            elif "THRESHOLD_PEAK_WIDTH_TAILING" == str_name:
                cfg.C5_THRESHOLD_PEAK_WIDTH_TAILING = float(str_value)

            elif "THRESHOLD_INVALID_ACQUIRING_SCAN" == str_name:
                cfg.C6_THRESHOLD_INVALID_ACQUIRING_SCAN = float(str_value)

            elif "FLAG_ANALYZE_FEATURE" == str_name:
                cfg.E3_FLAG_ANALYZE_FEATURE = int(str_value)

            elif "PATH_EXPORT" == str_name:      # 固定输出文件路径
                cfg.E1_PATH_EXPORT = str_value
                raw_file = cfg.A0_PATH_RAW.split('\\')[-1]
                raw_name = raw_file.split('.')[0]
                cfg.E1_PATH_EXPORT = cfg.E1_PATH_EXPORT + 'MSCohort_' + raw_name +'\\'
                if os.access(cfg.E1_PATH_EXPORT, os.F_OK):
                    pass
                else:
                    os.makedirs(cfg.E1_PATH_EXPORT)

            else:

                info = str_name + " in the config file is all right?"
                logGetWarning(info)


class CRaw2MS:

    def FunctionRaw2MS(self, inputListPathRaw):

        exe_MSFileReader = 'MSExport.exe'
        cfg_MSFileReader = 'MSExport.cfg'
        path_raw = '|'.join(inputListPathRaw)

        with open(cfg_MSFileReader, 'w') as f:

            f.write("[Data]\n")
            f.write("PATH_RAW=" + path_raw + '\n')

        cmd = exe_MSFileReader + ' {:s}'.format(cfg_MSFileReader)

        try:
            os.system(cmd)
        except:
            print('MSExport.exe run wrong!')


class CParseIDForSpectronaut:

    def __init__(self, inputDP):

        self.dp = inputDP

    def read(self, pathID):

        global position_PrecWindowNumber, position_organism
        with open(pathID, 'rb') as f:
            lines = f.read().decode('utf-8').split('\r\n')

        table_list = lines[0].split('\t')
        lines = lines[1:-1] if lines[-1] == '' else lines[1:]

        position_rawname = table_list.index('R.FileName')
        position_rundata = table_list.index('R.Run Date')
        position_gradient = table_list.index('R.Gradient Length [min]')
        position_PGQuant = table_list.index('PG.Quantity')
        position_PepQuant = table_list.index('PEP.Quantity')
        position_SignaltoNoise = table_list.index('EG.SignalToNoise')
        position_precursorQuant = table_list.index('EG.TotalQuantity (Settings)')


        position_proteingroup = table_list.index('PG.ProteinGroups')
        position_proteinqvalue = table_list.index('PG.Qvalue')
        position_precursorqvalue = table_list.index('EG.Qvalue')
        position_peptide = table_list.index('PEP.StrippedSequence')
        position_mod = table_list.index('EG.ModifiedPeptide')
        position_precursor = table_list.index('EG.PrecursorId')
        position_precursorRTStart = table_list.index('EG.StartRT')
        position_precursorRTEnd = table_list.index('EG.EndRT')
        position_precursorRT = table_list.index('EG.ApexRT')
        position_FWHM = table_list.index('EG.FWHM')
        position_PeakWidth = table_list.index('EG.PeakWidth')
        position_charge = table_list.index('FG.Charge')
        position_precursorMOZ = table_list.index('FG.PrecMz')
        position_MS1_MassAccuray = table_list.index('FG.RawMassAccuracy (PPM)')
        position_MS2_MassAccuray = table_list.index('F.RawMassAccuracy (PPM)')
        position_PrecWindow = table_list.index('FG.PrecWindow')
        position_PrecWindowNumber = table_list.index('FG.PrecWindowNumber')
        position_MissedCleavage = table_list.index('PEP.NrOfMissedCleavages')
        position_MS2_datapoint = table_list.index('EG.DatapointsPerPeak')
        position_MS1_datapoint = table_list.index('EG.DatapointsPerPeak (MS1)')
        position_decoy = table_list.index('EG.IsDecoy')
        last_ID = '  '

        for line in lines:

            line_list = line.split('\t')
            RAW_PATH = str(self.dp.myCFG.A0_PATH_RAW.split('\\')[-1])
            RAW_NAME = str(RAW_PATH.split('.')[0])
            tmp_rawname = line_list[position_rawname]
            tmp_decoy = line_list[position_decoy]

            if tmp_rawname == RAW_NAME and (tmp_decoy == 'FALSE' or tmp_decoy == 'False'):
                tmp_precursor = line_list[position_precursor]
                tmp_group = line_list[position_proteingroup]

                tmp_ID = tmp_rawname + '_' + tmp_group + '_' + tmp_precursor

                tmp_proteinScore = float(line_list[position_proteinqvalue])
                tmp_peptide = line_list[position_peptide]
                tmp_precursorMod = self.__soliderParseMod(line_list[position_mod])
                tmp_precursorCharge = int(line_list[position_charge])
                tmp_precursorRTStart = float(line_list[position_precursorRTStart])
                tmp_precursorRTEnd = float(line_list[position_precursorRTEnd])
                tmp_precursorRT = float(line_list[position_precursorRT])
                tmp_FWHM = float(line_list[position_FWHM])
                tmp_PeakWidth = float(line_list[position_PeakWidth])
                tmp_MS1_MassAccuray = float(line_list[position_MS1_MassAccuray])
                tmp_MS2_MassAccuray = float(line_list[position_MS2_MassAccuray])
                tmp_PrecWindow = line_list[position_PrecWindow]
                tmp_missedCleavage = int(line_list[position_MissedCleavage])
                tmp_MS2_datapoint = float(line_list[position_MS2_datapoint])
                tmp_MS1_datapoint = float(line_list[position_MS1_datapoint])
                tmp_PreWindowNumber = int(float(line_list[position_PrecWindowNumber]))

                tmp_precursorScore = float(line_list[position_precursorqvalue])
                tmp_precursorMOZ = float(line_list[position_precursorMOZ])
                tmp_PGQuant = float(line_list[position_PGQuant])

                tmp_PepQuant = float(line_list[position_PepQuant])
                tmp_precursorQuant = float(line_list[position_precursorQuant])
                tmp_SignaltoNoise = float(line_list[position_SignaltoNoise])

                if last_ID != tmp_ID:
                    last_ID = tmp_ID

                    self.dp.myID.PSM31_Precursor.append(tmp_precursor)
                    self.dp.myID.PSM12_SCORE2.append(tmp_precursorScore)
                    self.dp.myID.PSM15_PRE_MOZ.append(tmp_precursorMOZ)
                    self.dp.myID.PSM5_MOD.append(tmp_precursorMod)
                    self.dp.myID.PSM9_CHARGE.append(tmp_precursorCharge)
                    self.dp.myID.PSM3_RT.append(tmp_precursorRT)
                    self.dp.myID.PSM22_RT_START.append(tmp_precursorRTStart)
                    self.dp.myID.PSM23_RT_END.append(tmp_precursorRTEnd)
                    self.dp.myID.PSM28_FWHM.append(tmp_FWHM)
                    self.dp.myID.PSM29_PeakWidth.append(tmp_PeakWidth)
                    self.dp.myID.PSM32_DataPoint_MS1.append(tmp_MS1_datapoint)
                    self.dp.myID.PSM24_DataPoint.append(tmp_MS2_datapoint)
                    self.dp.myID.PSM25_PrecWindow.append(tmp_PrecWindow)
                    self.dp.myID.PSM25_PrecWindowNumber.append(tmp_PreWindowNumber)
                    self.dp.myID.PSM26_MS1_MassAccuracy.append(tmp_MS1_MassAccuray)
                    self.dp.myID.PSM27_MS2_MassAccuracy.append(tmp_MS2_MassAccuray)
                    self.dp.myID.PSM18_PRE_INTENSITY.append(tmp_precursorQuant)
                    self.dp.myID.PSM17_PRE_SignalToNoise.append(tmp_SignaltoNoise)
                    self.dp.myID.PSM30_MissedCleavage_PRE.append(tmp_missedCleavage)

                    self.dp.myID.N_PSM = self.dp.myID.N_PSM + 1

                    if tmp_peptide not in self.dp.myID.PSM4_SEQ:
                        self.dp.myID.PSM4_SEQ.append(tmp_peptide)
                        self.dp.myID.PSM14_PEP_INTENSITY.append(tmp_PepQuant)
                        self.dp.myID.PSM30_MissedCleavage.append(tmp_missedCleavage)

                    if tmp_group not in self.dp.myID.PSM20_PRO:
                        self.dp.myID.PSM20_PRO.append(tmp_group)
                        self.dp.myID.PSM11_SCORE1.append(tmp_proteinScore)
                        self.dp.myID.PSM13_PRO_INTENSITY.append(tmp_PGQuant)

                    if tmp_rawname not in self.dp.myID.PSM1_RAW_NAME:
                        self.dp.myID.PSM1_RAW_NAME.append(tmp_rawname)
                        self.dp.myID.R1_RUN_DATA.append(line_list[position_rundata])
                        self.dp.myID.R2_GRADIENT.append(float(line_list[position_gradient]))


                    if 'CON' not in tmp_group and 'REV' not in tmp_group:
                        if ';' in tmp_group:
                            tmp_group = tmp_group.split(';')
                            for i_pro in tmp_group:
                                self.dp.myID.PSM20_PROTEIN.append(i_pro)
                        else:
                            self.dp.myID.PSM20_PROTEIN.append(tmp_group)
        self.dp.myID.PSM20_PROTEIN = set(self.dp.myID.PSM20_PROTEIN)


    def __soliderParseMod(self, inputMod):

        str_mod = ''
        flag_mod = 0

        for i in inputMod:
            if i == '[':
                flag_mod = 1
            elif i == ']':
                flag_mod = 0
                str_mod += ';'
            elif flag_mod == 1:
                str_mod += i

        return str_mod


class CAnalyzeIDMassDeviation:

    def __init__(self, inputDP):
        self.dp = inputDP

    def analyze(self):
        path = self.dp.myCFG.E1_PATH_EXPORT + '\\' + IO_FILENAME_EXPORT[4]
        if self.dp.myCFG.B3_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['Spectronaut']:
            with open(path, 'w') as f:
                f.write('\t'.join(["PrecursorID", "RT", "Precursor_Qvalue", "MS1_MassAccuracy", "MS2_MassAccuracy"]) + '\n')
                for i in range(self.dp.myID.N_PSM):
                    f.write(str(self.dp.myID.PSM31_Precursor[i]))
                    f.write('\t')
                    f.write('{0:.5f}'.format(self.dp.myID.PSM3_RT[i]))
                    f.write('\t')
                    f.write('{0:e}'.format(self.dp.myID.PSM12_SCORE2[i]))
                    f.write('\t')
                    f.write('{0:.5f}'.format(self.dp.myID.PSM26_MS1_MassAccuracy[i]))
                    f.write('\t')
                    f.write('{0:.5f}'.format(self.dp.myID.PSM27_MS2_MassAccuracy[i]))
                    f.write('\n')


class CAnalyzeDIAchromatography:

    def __init__(self, inputDP):
        self.dp = inputDP

    def analyze(self):

        path = self.dp.myCFG.E1_PATH_EXPORT + '\\' + IO_FILENAME_EXPORT[5]
        with open(path, 'w') as f:
            f.write('\t'.join(
                ["PrecursorID", "RT_START", "RT_END", "Peak_Width", "FWHM", "DataPoint_MS2", "DataPoint_MS1"]) + '\n')
            n_psm = self.dp.myID.N_PSM

            for i_psm in range(n_psm):
                f.write(str(self.dp.myID.PSM31_Precursor[i_psm]))
                f.write('\t')
                f.write('{0:.2f}'.format(self.dp.myID.PSM22_RT_START[i_psm]))
                f.write('\t')
                f.write('{0:.2f}'.format(self.dp.myID.PSM23_RT_END[i_psm]))
                f.write('\t')
                f.write('{0:.3f}'.format(self.dp.myID.PSM29_PeakWidth[i_psm]))
                f.write('\t')
                f.write('{0:.3f}'.format(self.dp.myID.PSM28_FWHM[i_psm]))
                f.write('\t')
                f.write('{0:.2f}'.format(self.dp.myID.PSM24_DataPoint[i_psm]))
                f.write('\t')
                f.write('{0:.2f}'.format(self.dp.myID.PSM32_DataPoint_MS1[i_psm]))
                f.write('\n')

        path2 = self.dp.myCFG.E1_PATH_EXPORT + '\\' + IO_FILENAME_EXPORT[3]
        with open(path2, 'w') as f2:
            f2.write('\t'.join(
                ["PrecursorID", "RT", "Precursor_Mz", "EG_TargetQuantity", "EG_SignalToNoise", "PrecWindow", "PrecWindowNumber","NrOfMissedCleavages"]) + '\n')
            n_psm = self.dp.myID.N_PSM

            for i_psm in range(n_psm):
                f2.write(str(self.dp.myID.PSM31_Precursor[i_psm]))
                f2.write('\t')
                f2.write('{0:.2f}'.format(self.dp.myID.PSM3_RT[i_psm]))
                f2.write('\t')
                f2.write('{0:.2f}'.format(self.dp.myID.PSM15_PRE_MOZ[i_psm]))
                f2.write('\t')
                f2.write('{0:.2f}'.format(self.dp.myID.PSM18_PRE_INTENSITY[i_psm]))
                f2.write('\t')
                f2.write('{0:.2f}'.format(self.dp.myID.PSM17_PRE_SignalToNoise[i_psm]))
                f2.write('\t')
                f2.write(str(self.dp.myID.PSM25_PrecWindow[i_psm]))
                f2.write('\t')
                f2.write(str(self.dp.myID.PSM25_PrecWindowNumber[i_psm]+1))
                f2.write('\t')
                f2.write(str(self.dp.myID.PSM30_MissedCleavage_PRE[i_psm]))
                f2.write('\n')
















#!/usr/bin/env python
# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import shutil, os, time, sys, filecmp
import dicom
#count the directories in a folder
def countDirectories(folderName):
    subDirectories = os.listdir(folderName)
    return len(subDirectories)
#read the RP file for the ptid, ptLastName,ptFirstName,planName
def readRPfile(fileName):
    dataSet=dicom.read_file(fileName)
    patientID=dataSet.PatientID
    patientName=dataSet.PatientName
    planName=dataSet.RTPlanLabel
    #this may not be stable may need to use the OS creation date and time
    planDateTime=dataSet.InstanceCreationDate+"-"+dataSet.InstanceCreationTime
    patientLastName=patientName.split('^')[0]
    patientFirstName=patientName.split('^')[1]
    return patientID, patientLastName, patientFirstName, planName, planDateTime
#Move data to processed folder
def createDirName(folderName,ptId,ptLastName,ptFirstName,ptPlanName,ptPlanDateTime):
    subDirectory="-"
    subDirectory=subDirectory.join([ptId,ptLastName,ptFirstName,ptPlanName,ptPlanDateTime])
    return subDirectory
#Move data to processed folder continued
def createDir(folderName,subDirectory):
    directory=folderName+subDirectory
    if not os.path.exists(directory):
        os.makedirs(directory)
        return 0
    else:
        return 1
        
def moveFiles(filesIn,sourceName,destName):
    for files in filesIn:
        shutil.move(sourceName+files,destName)
#walk through the directories and process each one
def processDirectories(folderInName,folderOutName):
    errorCode=[]
    skipDir=[]
    subDirectories = os.listdir(folderInName)
    #check is there are CT, RP, RS, RD files
    for subDirs in subDirectories:
        #print subDirs
        files=os.listdir(folderInName+subDirs)
        CTfiles=[];RSfiles=[];RPfiles=[];RDfiles=[];
        for ifile in files:
            spiFile=ifile.split(".")[0]
            #print ifile
            if spiFile=="CT":
                CTfiles.append(ifile)
            elif spiFile=="RS":
                RSfiles.append(ifile)
            elif spiFile=="RP":
                RPfiles.append(ifile)
            elif spiFile=="RD":
                RDfiles.append(ifile)
        errorCode.append(subDirs+" ")
        skipDir.append(0)
        errorCodeLen=len(errorCode)
        index=errorCodeLen-1
        if not len(CTfiles):
            errorCode[index]=errorCode[index]+"No CT files: "
        if not len(RSfiles):
            errorCode[index]=errorCode[index]+"No RS files: "
        if not len(RPfiles):
            errorCode[index]=errorCode[index]+"No RP files: "
        if not len(RDfiles):
            errorCode[index]=errorCode[index]+"No RD files: "
        #add check for newest files
        #update them so that only newest files are used
        if not errorCode[index]==subDirs+" ":
            #print errorCode[index]
            #print 'Skipping directory '+subDirs
            skipDir[index]=1
        if not skipDir[index]:
            #print 'Processing directory '+subDirs
            ptId,ptLastName,ptFirstName,ptPlanName,ptPlanDateTime=readRPfile(folderInName+subDirs+"/"+RPfiles[0])
            #print ptId, ptLastName, ptFirstName, ptPlanName
            subDirectory=createDirName(folderOutName,ptId,ptLastName,ptFirstName,ptPlanName,ptPlanDateTime)
            #print subDirectory
            if not createDir(folderOutName,subDirectory):
                moveFiles(CTfiles,folderInName+subDirs+"/",folderOutName+subDirectory)
                moveFiles(RSfiles,folderInName+subDirs+"/",folderOutName+subDirectory)
                moveFiles(RPfiles,folderInName+subDirs+"/",folderOutName+subDirectory)
                moveFiles(RDfiles,folderInName+subDirs+"/",folderOutName+subDirectory)
                #remove folderInName
                shutil.rmtree(folderInName+subDirs)
                errorCode[index]="Successfully processed "+subDirs
            else:
                #compare the files since the directory already existed
                dcmp=filecmp.dircmp(folderInName+subDirs,folderOutName+subDirectory)
                print 'Length of diff_files ',len(dcmp.left_only)
                if len(dcmp.left_only)==0:
                    #all files are the same just delete the ones in DICOMS and move on
                    shutil.rmtree(folderInName+subDirs)
                    errorCode[index]="Successfully processed "+subDirs
                else:
                    #copy the different files
                    print 'Updating new files'
                    for name in dcmp.left_only:
                        shutil.move(folderInName+subDirs+"/"+name,folderOutName+subDirectory)
                    shutil.rmtree(folderInName+subDirs)
                    errorCode[index]="Successfully processed "+subDirs
    return errorCode,skipDir
    
DICOMSdir='/home/emcuser/eMCfiles/DICOMS/'
DICOMSProcdir='/home/emcuser/eMCfiles/DICOMSPROC/'
RUNNER=True
def main(RUNNER):
    k=0    
    while RUNNER:
        try:
            print countDirectories(DICOMSdir)
            countDirs=countDirectories(DICOMSdir)
            #if more than zero search each one for patient id, lastname, firstname, planname
            if countDirs > 0:
                errors,skips=processDirectories(DICOMSdir, DICOMSProcdir)
                if not len(errors)==0:
                    for error in errors:
                        print error
        except Exception,e:
            print str(e)
            print 'stopping'
            #RUNNER=False
            log_file=open('/home/emcuser/eMCfiles/log/process.log','a')
            log_file.write(str(e)+"\n")
            log_file.close()
        k=k+1
        print 'current iteration ',k
        time.sleep(120)
    sys.exit(0)
if __name__ == "__main__":
    main(RUNNER)
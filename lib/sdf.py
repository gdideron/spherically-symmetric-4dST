# Module for reading in sdf files 
# W. East, 2013

import struct
from collections import namedtuple
import ctypes as ct
import numpy as np
import sys


# Error classes    
class Error(Exception):
    pass

class HeaderError(Error):
    pass

class DataCorruptError(Error):
    pass

class TimeStepInconsitenceError(Error):
    pass


class SDFData(dict):
    """
    SDF file for manipulation.

    Attributes:
        time             : Initial time
        timeArray        : All data times
        rank             : Dimension of data
        dataName
        data
        coordName
        coordData
    """

    def __init__(self, fileStream, loadData=True):
        header = self._loadHeader(fileStream)
        self.time = header.time
        self.timeArray = np.array([])
        self.rank = header.rank
        self._version = header.version
        self._csize = header.csize
        self._dsize = header.dsize
        self._pnlen = header.pnlen
        self._cnlen = header.cnlen
        self._tglen= header.tglen

        self._loadVars(fileStream)

        if (loadData):
            while True:
                self.timeArray = np.hstack((self.timeArray,header.time))
                self._loadData(fileStream)
                try:
                    header = self._loadHeader(fileStream)
                except HeaderError:
                    break
                self._loadVars(fileStream)
        else:
            #Save information to load data later 
            self._fileName = fileStream.name
            self._fileTell = fileStream.tell() 
            # Seek past data
            fileStream.seek(ct.sizeof(ct.c_double)*self._csize, 1)
            fileStream.seek(ct.sizeof(ct.c_double)*self._dsize, 1)
            self._data = np.array([])
            self._coordData = np.array([])

    def _loadHeader(self, fileStream):
        """
        Giving a fileStream pointing to the start of a rnpletal header, reads
        the header into the object attributes.

        """
        fmt = '!dddddddd' # Big-Endian 8 c-type doubles
        head = fileStream.read(struct.calcsize(fmt))
        if not head:
            raise HeaderError
        FileHeader = namedtuple("FileHeader", 
                            "time version rank dsize csize pnlen cnlen tglen")
        try:
            ( time, version, rank, dsize, 
              csize, pnlen, cnlen, tglen ) = struct.unpack(fmt,head)
        except struct.error: 
            raise HeaderError
        header = FileHeader( time, int(version), int(rank), int(dsize),
                             int(csize), int(pnlen), int(cnlen), int(tglen) )
        return header

    def _loadVars(self,fileStream):
        """ Loads data name, coordinate name, boundary box and shape"""
        try:
            dataName = fileStream.read(self._pnlen*ct.sizeof(ct.c_char))
        except OverflowError:
            raise HeaderError
        try:
            if isinstance(self.dataName,str):
                if (dataName != self.dataName):
                    raise TimeStepInconsitenceError
        except AttributeError:
            self.dataName = dataName
            if (b"_c_" in self.dataName):
                self.isCellCentered = True
            else: 
                self.isCellCentered = False

        try:
            coordName = fileStream.read(self._cnlen*ct.sizeof(ct.c_char))
        except OverflowError:
            raise HeaderError 
        try:
            if isinstance(self.coordName,str):
                if (coordName!= self.coordName):
                    raise TimeStepInconsitenceError
        except AttributeError:
            self.coordName = coordName 

        dt = np.dtype(np.float64).newbyteorder('B')
        try:
            bbox = np.fromfile(fileStream, dt, 2*self.rank)
            shape = np.fromfile(fileStream, dt, self.rank).astype(int)
        except OverflowError:
            raise HeaderError 

        try:
            if isinstance(self.bbox,np.ndarray):
                np.testing.assert_almost_equal(self.bbox,bbox)
            if isinstance(self.shape,np.ndarray):
                np.testing.assert_almost_equal(self.shape,shape)
        except AttributeError:
           self.bbox = bbox 
           self.shape = shape 
           range1 = range(0,2*self.rank,2)
           range2 = range(1,2*self.rank,2)
           if (self.rank*2!=self.bbox.size or self.rank!=self.shape.size): 
               raise HeaderError
           self._dx = (self.bbox[range2] - self.bbox[range1])/(self.shape-1.0)
           if (self.isCellCentered): 
               #Adjust for the way PAMR saves the bbox as ending at cell centers
               self.bbox[range1] = self.bbox[range1] - 0.5*self._dx
               self.bbox[range2] = self.bbox[range2] + 0.5*self._dx

    def _loadData(self, fileStream):
        dt = np.dtype(np.float64).newbyteorder('B')
        if (self._csize>0):
            try:
                coordData = np.fromfile(fileStream, dt, self._csize)
            except OverflowError:
                print("Data is corrupted")
                raise DataCorruptError
                data = 0 
        else:
            coordData = np.array([])
        try:
            if isinstance(self._coordData,np.ndarray):
                np.testing.assert_almost_equal(self._coordData,coordData)
        except AttributeError:
           self._coordData = coordData 

        try:
            data = np.fromfile(fileStream, dt, self._dsize)
        except OverflowError:
            print("Data is corrupted")
            raise DataCorruptError
            data = 0 
        if (data.size!=np.prod(self.shape)):
            print("Data is corrupted")
            raise DataCorruptError
            data = 0 
        else:
            data = data.reshape(self.shape[::-1]).T
        try:
            if isinstance(self._data,np.ndarray):
                self._data = np.vstack((self._data,data))
        except AttributeError:
            self._data = np.array([data])

    def getData(self):
        if (not self._data.size):
            with open(self._fileName, "rb") as f:
                f.seek(self._fileTell)
                self._loadData(f)
                f.close()
        return self._data 

    def getDataOnce(self):
        with open(self._fileName, "rb") as f:
            f.seek(self._fileTell)
            dt = np.dtype(np.float64).newbyteorder('B')
            if (self._csize>0):
                coordData = np.fromfile(f, dt, self._csize)
            else:
                coordData = np.array([])
            dataf = np.fromfile(f, dt, self._dsize)
            if (dataf.size!=np.prod(self.shape)):
                print("Data is corrupted")
                raise DataCorruptError
                data = 0
            else:
                data = dataf.reshape(self.shape[::-1]).T
            f.close()
        return data

    def setData(self, data):
        if (self._data.size==data.size):
            dt = np.dtype(np.float64).newbyteorder('B')
            self._data = data.astype(dt) 
        else:
            print("Error in setData: data size doesn't match.")

    def freeData(self):
        if (not self._data.size):
            del self._data
            self._data = np.array([])
        if (not self._coordData.size):
            del self._coordData
            self._coordData = np.array([])

    def writeToFile(self, fileStream):
        fmt = '!dddddddd'
        self._data = self.getData()
        for i in range(self.timeArray.size):
            time = self.timeArray[i]
            head = struct.pack( fmt, time, float(self._version), 
                                float(self.rank), float(self._dsize), 
                                float(self._csize), float(self._pnlen),
                                float(self._cnlen), float(self._tglen) )
            fileStream.write(head)
            fileStream.write(self.dataName)
            fileStream.write(self.coordName)
            self.bbox.tofile(fileStream)
            dt = np.dtype(np.float64).newbyteorder('B')
            self.shape.astype(dt).tofile(fileStream)
            data = (self._data[i].T).reshape(self.shape[::-1]).astype(dt)
            if (self._csize>0):
                self._coordData.tofile(fileStream)
            data.tofile(fileStream)

    def getDx(self):
        if self.rank<1:
            return 0
        return self._dx[0]      

    def getCoord(self, rank=0):
        if (self._coordData.size == np.sum(self.shape)):
            i = 0
            for n in range(rank):
                i = i + self.shape[n]
            return self._coordData[i:i+self.shape[rank]]
        return np.linspace(self.bbox[2*rank], self.bbox[2*rank+1], num=self.shape[rank])

# Returns a list of lists of sdfData objects at the same time
def groupByTime(sdfFiles):
    sdfFiles = sortByTime(sdfFiles)
    times = []
    sdfFilesAtTimes = []
    for f in sdfFiles:
        if (times and np.allclose(f.time, times[-1])):
            (sdfFilesAtTimes[-1]).append(f)
        else:
            sdfFilesAtTimes.append([f])
            times.append(f.time)
    return [sdfFilesAtTimes, times]

def loadSDFFiles(fileNames, loadData=True):
    sdfFiles = []
    for fileName in fileNames:
        with open(fileName, "rb") as f:
            i=0
            while True:
                try:
                    sdfData = SDFData(f, loadData)
                except (HeaderError,DataCorruptError): break
                sdfFiles.append(sdfData)
        f.close()
    return sdfFiles
 
def purgeBadFiles(sdfFiles):
    purgedSdfFiles = []
    for f in sdfFiles:
        try:
            d = f.getData()
        except (HeaderError,DataCorruptError): continue 
        if (np.isnan(np.min(d))==False and np.max(abs(d))<1.0e20):
            purgedSdfFiles.append(f)
    return purgedSdfFiles

def writeSDFToFiles(sdfFiles, fileName): 
    with open(fileName, "wb") as f:
        for sdf in sdfFiles:
            sdf.writeToFile(f)

#Takes a numpy array and writes sdf to filestream
def writeArrayToSDF( data, bbox, fileStream, dataName="data", 
                     coordName="x", coordData=None, timeArray=0.0 ):
    try:
        len(timeArray)
    except TypeError:
        timeArray = [timeArray]
    if len(data.shape) == 1:
        dataStruct = data
        data = [data]
    else:
        dataStruct = data[0]
    #Need to make strings null terminating
    if (len(dataName)==0 or dataName[-1]!='\0'):
        dataName = dataName+'\0'
    if (len(coordName)==0 or coordName[-1]!='\0'):
        coordName = coordName+'\0'
    if (coordData is None):
        coordData = bbox #Seems to cause problems if coordData is empty
    dsize = dataStruct.size 
    shape = np.array(dataStruct.shape)
    rank = len(shape)
    csize = coordData.size
    pnlen = len(dataName)
    cnlen = len(coordName)
    tglen = 0 
    version = 1
    fmt = '!dddddddd'
    for i, time in enumerate(timeArray):
        head = struct.pack( fmt, time, float(version), float(rank),
                            float(dsize), float(csize), float(pnlen),
                            float(cnlen), float(tglen) )
        fileStream.write(head)
        if sys.version_info < (3, 0):
            fileStream.write(dataName)
            fileStream.write(coordName)
        else:
            fileStream.write(bytes(dataName,'utf-8'))
            fileStream.write(bytes(coordName,'utf-8'))
        dt = np.dtype(np.float64).newbyteorder('B')
        bbox.astype(dt).tofile(fileStream)
        shape.astype(dt).tofile(fileStream)
        dataIn = (data[i].T).reshape(shape[::-1])
        coordData.astype(dt).tofile(fileStream)
        dataIn.astype(dt).tofile(fileStream)
